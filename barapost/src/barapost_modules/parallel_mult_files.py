# -*- coding: utf-8  -*-
# This module defines functions necessary for barapost.py to perform parallel
#   processing in "many-files" mode.

import os
import shelve
import multiprocessing as mp
from re import search as re_search

from src.fasta import fasta_packets
from src.fastq import fastq_packets

from src.printlog import getwt, printl, printn, println
from src.write_classification import write_classification
from src.spread_files_equally import spread_files_equally
from src.filesystem import get_curr_res_dpath, create_result_directory
from src.filesystem import remove_tmp_files, is_fastq, OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped

from src.barapost_modules.barapost_spec import look_around, launch_blastn, parse_align_results_xml


def init_process(print_lock_buff, conter_lock_buff, file_counter_buff):
    """
    Function initializes global variables that all processes shoud have access to.
    This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.

    :param print_lock_buff: lock that synchronizes printing to the console;
    :type print_lock_buff: multiprocessing.Lock;
    :param counter_lock_buff: lock that synchronizes incrementing 'file_counter';
    :type counter_lock_buff: multiprocessing.Lock;
    :param file_counter: variable for counting processed files;
    :type file_cunter: mp.Value('i');
    """

    global print_lock
    print_lock = print_lock_buff

    global counter_lock
    counter_lock = conter_lock_buff

    global file_counter
    file_counter = file_counter_buff
# end def init_proc_many_files


def process_paral(fq_fa_list, packet_size, tax_annot_res_dir,
    blast_algorithm, use_index, db_path, nfiles, logfile_path):
    """
    Function performs 'many_files'-parallel mode of barapost.py.

    :param fq_fa_list: list of paths to FASTA and FASTQ files meant to be processed;
    :type fq_fa_list: list<str>;
    :param packet_size: number of sequences processed by blast in a single launching;
    :type packet_size: int;
    :param tax_annot_res_dir: path to ouput directory that contains taxonomic annotation;
    :type tax_annot_res_dir: str;
    :param blast_algorithm: blast algorithm to use;
    :type blast_algorithm: str;
    :param use_index: logic value indicationg whether to use indes;
    :type use_index: bool;
    :param db_path: path to database;
    :type db_path: str;
    :param nfiles: total number of files;
    :type nfiles: int;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    queries_tmp_dir = os.path.join(tax_annot_res_dir, "queries-tmp")

    # Iterate over source FASTQ and FASTA files
    for fq_fa_path in fq_fa_list:

        # Create the result directory with the name of FASTQ of FASTA file being processed:
        new_dpath = create_result_directory(fq_fa_path, tax_annot_res_dir)

        # "hname" means human readable name (i.e. without file path and extention)
        infile_hname = os.path.basename(fq_fa_path)
        infile_hname = re_search(r"(.+)\.(m)?f(ast)?(a|q)(\.gz)?$", infile_hname).group(1)

        # Look around and ckeck if there are results of previous runs of this script
        # If 'look_around' is None -- there is no data from previous run
        previous_data = look_around(new_dpath, fq_fa_path, blast_algorithm, logfile_path)

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
            num_seqs = len(tuple(filter(lambda l: True if l.startswith('>') else False,
                map(fmt_func, how_to_open(fq_fa_path).readlines()))))
        # end if

        if num_seqs == num_done_seqs:
            with counter_lock:
                file_counter.value += 1
                i = file_counter.value # save to local var and release lock
            # end with
            with print_lock:
                printl(logfile_path, "\r{} - File #{}/{} ('{}') has been already completely processed.".format(getwt(),
                    i, nfiles, fq_fa_path))
                println(logfile_path, "Omitting it.\nWorking...")
            # end with
            continue
        # end if

        for packet in packet_generator(fq_fa_path, packet_size, num_done_seqs):

            # Blast the packet
            align_xml_text = launch_blastn(packet["fasta"], blast_algorithm,
                use_index, queries_tmp_dir, db_path, logfile_path)

            # Cnfigure result TSV lines
            result_tsv_lines = parse_align_results_xml(align_xml_text,
                packet["qual"], logfile_path)

            # Write the result to tsv
            write_classification(result_tsv_lines, tsv_res_path)
        # end for

        with counter_lock:
            file_counter.value += 1
            i = file_counter.value # save to local var and release lock
        # end with
        with print_lock:
            println(logfile_path, "\r{} - File #{}/{} ({}) is processed.\nWorking...".format(getwt(),
                i, nfiles, os.path.basename(fq_fa_path)))
        # end with
    # end for

    remove_tmp_files( os.path.join(queries_tmp_dir, "query{}_tmp.fasta".format(os.getpid())) )
# end def process_paral


def process(fq_fa_list, n_thr, packet_size, tax_annot_res_dir,
    blast_algorithm, use_index, db_path, logfile_path):
    """
    Function launches parallel processing in "many-files" mode by barapost.py.

    :param fq_fa_list: list of paths to files meant to be processed;
    :type fq_fa_list: list<str>;
    :param n_thr: number of threads to launch;
    :type n_thr: int;
    :param packet_size: number of sequences processed by blast in a single launching;
    :type packet_size: int;
    :param tax_annot_res_dir: path to ouput directory that contains taxonomic annotation;
    :type tax_annot_res_dir: str;
    :param blast_algorithm: blast algorithm to use;
    :type blast_algorithm: str;
    :param use_index: logic value indicationg whether to use indes;
    :type use_index: bool;
    :param db_path: path to database;
    :type db_path: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    pool = mp.Pool(n_thr, initializer=init_process,
        initargs=(mp.Lock(), mp.Lock(), mp.Value('i', 0),))

    pool.starmap(process_paral, [ (fq_fa_sublist,
        packet_size,
        tax_annot_res_dir,
        blast_algorithm,
        use_index,
        db_path,
        len(fq_fa_list),
        logfile_path) for fq_fa_sublist in spread_files_equally(fq_fa_list, n_thr) ])

    # Reaping zombies
    pool.close()
    pool.join()
# end def process
