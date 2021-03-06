# -*- coding: utf-8 -*-
# This module defines functions necessary for barapost-local.py to perform parallel
#   processing in "few-files" mode.

import os
import re
import sys
import multiprocessing as mp

from src.fasta import fasta_packets, fasta_packets_from_str
from src.fastq import fastq_packets

from src.printlog import printlog_info, printlog_info_time, printn, printlog_warning
from src.write_classification import write_classification
from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped, is_fastq
from src.filesystem import create_result_directory, remove_tmp_files

from src.barapost_local_modules.barapost_spec import look_around, launch_blastn, parse_align_results_xml


def init_proc_single_file_in_paral(print_lock_buff, write_lock_buff):
    # Function initializes global variables that all processes shoud have access to.
    # This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.

    # :param print_lock_buff: lock that synchronizes printing to the console;
    # :type print_lock_buff: multiprocessing.Lock;
    # :param write_lock_buff: lock that synchronizes wriiting to result file;
    # :type write_lock_buff: multiprocessing.Lock;

    global print_lock
    print_lock = print_lock_buff

    global write_lock
    write_lock = write_lock_buff
# end def init_proc_single_file_in_paral


def process_part_of_file(data, tsv_res_path, packet_size, tax_annot_res_dir,
    blast_algorithm, use_index, db_path):
    # Function preforms processing part of file in 'few_files'-parallel mode.

    # :param data: fasta-formatted string meant to be processed;
    # :type data: str;
    # :param tsv_res_path: path to result TSV file;
    # :type tsv_res_path: str;
    # :param packet_size: number of sequences processed by blast in a single launching;
    # :type packet_size: int;
    # :param tax_annot_res_dir: path to ouput directory that contains taxonomic annotation;
    # :type tax_annot_res_dir: str;
    # :param blast_algorithm: blast algorithm to use;
    # :type blast_algorithm: str;
    # :param use_index: logic value indicationg whether to use indes;
    # :type use_index: bool;
    # :param db_path: path to database;
    # :type db_path: str;

    queries_tmp_dir = os.path.join(tax_annot_res_dir, "queries-tmp")

    for packet in fasta_packets_from_str(data["fasta"], packet_size):

        # Blast the packet
        align_xml_text = launch_blastn(packet["fasta"], blast_algorithm,
            use_index, queries_tmp_dir, db_path)

        # Get result tsv lines
        result_tsv_lines = parse_align_results_xml(align_xml_text, data["qual"])
        # If we use packet["qual"] -- we will have all '-'-s because 'data' is a fasta-formatted string
        # Thus there are no value for key "qual" in 'packet' (see src/barapost_local_modules/fasta_packets_from_str.py)

        # Write the result to TSV file
        with write_lock:
            write_classification(result_tsv_lines, tsv_res_path)
        # end with
    # end for

    query_fpath = os.path.join(queries_tmp_dir, "query{}_tmp.fasta".format(os.getpid()))
    remove_tmp_files(query_fpath)
# end def process_part_of_file


def process(fq_fa_list, n_thr, packet_size, tax_annot_res_dir,
            blast_algorithm, use_index, db_path):
    # Function preforms "few_files"-parallel mode.
    #
    # :param fq_fa_list: list of paths to files meant to be processed;
    # :type fq_fa_list: list<str>;
    # :param n_thr: number of threads to launch;
    # :type n_thr: int;
    # :param packet_size: number of sequences processed by blast in a single launching;
    # :type packet_size: int;
    # :param tax_annot_res_dir: path to ouput directory that contains taxonomic annotation;
    # :type tax_annot_res_dir: str;
    # :param blast_algorithm: blast algorithm to use;
    # :type blast_algorithm: str;
    # :param use_index: logic value indicationg whether to use indes;
    # :type use_index: bool;
    # :param db_path: path to database;
    # :type db_path: str;

    nfiles = len(fq_fa_list)

    for i, fq_fa_path in enumerate(fq_fa_list):
        # Create the result directory with the name of FASTQ of FASTA file being processed:
        new_dpath = create_result_directory(fq_fa_path, tax_annot_res_dir)

        # "hname" means human readable name (i.e. without file path and extention)
        infile_hname = os.path.basename(fq_fa_path)
        infile_hname = re.search(r"(.+)\.(m)?f(ast)?(a|q)(\.gz)?$", infile_hname).group(1)

        # Look around and ckeck if there are results of previous runs of this script
        # If 'look_around' is None -- there is no data from previous run
        previous_data = look_around(new_dpath, fq_fa_path)

        if previous_data is None: # If there is no data from previous run
            num_done_seqs = 0 # number of successfully processed sequences
            tsv_res_path = os.path.join(new_dpath, "classification.tsv") # form result tsv file path
        else: # if there is data from previous run
            num_done_seqs = previous_data["n_done_reads"] # get number of successfully processed sequences
            tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        # end if

        how_to_open = OPEN_FUNCS[ is_gzipped(fq_fa_path) ]
        fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_fa_path) ]

        if is_fastq(fq_fa_path):
            packet_generator = fastq_packets
            num_seqs = sum(1 for line in how_to_open(fq_fa_path)) // 4 # 4 lines per record
        else:
            packet_generator = fasta_packets
            try:
                num_seqs = len(tuple(filter(lambda l: True if l.startswith('>') else False,
                    map(fmt_func, how_to_open(fq_fa_path).readlines()))))
            except UnicodeDecodeError as err:
                print()
                printlog_warning("Warning: current file is broken: {}."\
                    .format(str(err)))
                printlog_warning("File: `{}`".format(os.path.abspath(fq_fa_path)))
                printlog_warning("This file will not be processed.")
                continue
            # end try
        # end if

        packet_size = min(packet_size, num_seqs // n_thr)

        if num_seqs == num_done_seqs:
            sys.stdout.write('\r')
            printlog_info_time("File #{}/{} (`{}`) has been already completely processed."\
                .format(i+1, nfiles, fq_fa_path))
            printlog_info("Omitting it.")
            printn("Working...")
            return
        # end if

        # Get number of seqeunces to pass to each thread
        file_part_size = num_seqs // n_thr
        if num_seqs % n_thr != 0:
            file_part_size += 1
        # end if

        pool = mp.Pool(n_thr, initializer=init_proc_single_file_in_paral,
            initargs=(mp.Lock(), mp.Lock(),))

        pool.starmap(process_part_of_file, [(file_part,
            tsv_res_path,
            packet_size,
            tax_annot_res_dir,
            blast_algorithm,
            use_index,
            db_path) for file_part in packet_generator(fq_fa_path,
                file_part_size,
                num_done_seqs)])

        # Reaping zombies
        pool.close()
        pool.join()

        sys.stdout.write('\r')
        printlog_info_time("File #{}/{} (`{}`) is processed.".\
            format(i+1, nfiles, os.path.basename(fq_fa_path)))
        printn("Working...")    # end for
# end def process
