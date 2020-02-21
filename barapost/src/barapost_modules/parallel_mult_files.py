# -*- coding: utf-8  -*-
# This module defines functions necessary for barapost.py to perform parallel
#   processing in "many-files" mode.

import os
import multiprocessing as mp
from re import search as re_search

from src.fasta import fasta_packets
from src.fastq import fastq_packets

from src.printlog import getwt, printl, printn
from src.write_classification import write_classification
from src.spread_files_equally import spread_files_equally
from src.filesystem import get_curr_res_dpath, create_result_directory
from src.filesystem import remove_tmp_files, is_fastq, OPEN_FUNCS, is_gzipped

from src.barapost_modules.barapost_spec import look_around, launch_blastn, parse_align_results_xml


def init_process(print_lock_buff, packet_size_buff, tax_annot_res_dir_buff,
    blast_algorithm_buff, use_index_buff, logfile_path_buff):
    """
    Function initializes global variables that all processes shoud have access to.
    This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.

    :param print_lock_buff: lock that synchronizes printing to the console;
    :type print_lock_buff: multiprocessing.Lock;
    :param packet_size_buff: number of sequences processed by blast in a single launching;
    :type packet_size_buff: int;
    :param tax_annot_res_dir_buff: path to ouput directory that contains taxonomic annotation;
    :type tax_annot_res_dir_buff: str;
    :param blast_algorithm_buff: blast algorithm to use;
    :type blast_algorithm_buff: str;
    :param use_index_buff: logic value indicationg whether to use indes;
    :type use_index_buff: bool;
    :param logfile_path_buff: path to log file;
    :type logfile_path_buff: str;
    """

    global print_lock
    print_lock = print_lock_buff

    global packet_size
    packet_size = packet_size_buff

    global tax_annot_res_dir
    tax_annot_res_dir = tax_annot_res_dir_buff

    global blast_algorithm
    blast_algorithm = blast_algorithm_buff

    global use_index
    use_index = use_index_buff

    global logfile_path
    logfile_path = logfile_path_buff
# end def init_proc_many_files


def process_paral(fq_fa_list):
    """
    Function performs 'many_files'-parallel mode of barapost.py.

    :param fq_fa_list: list of paths to FASTA and FASTQ files meant to be processed;
    :type fq_fa_list: list<str>;
    """

    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy","taxonomy")
    queries_tmp_dir = os.path.join(tax_annot_res_dir, "queries-tmp")
    local_fasta = os.path.join(tax_annot_res_dir, "local_database", "local_seq_set.fasta")

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

        if is_fastq(fq_fa_path):
            packet_generator = fastq_packets
            num_seqs = sum(1 for line in how_to_open(fq_fa_path)) // 4 # 4 lines per record
        else:
            packet_generator = fasta_packets
            num_seqs = len(tuple(filter(lambda l: True if l.startswith('>') else False,
                map(fmt_func, how_to_open(fq_fa_path).readlines()))))
        # end if

        if num_seqs == num_done_seqs:
            with print_lock:
                printl(logfile_path, "\rFile '{}' has been already completely processed.".format(fq_fa_path))
                printl(logfile_path, "Omitting it.")
                printn("  Working...")
            # end with
            continue
        # end if

        for packet in packet_generator(fq_fa_path, packet_size, num_done_seqs):

            # Blast the packet
            align_xml_text = launch_blastn(packet["fasta"], blast_algorithm,
                use_index, queries_tmp_dir, local_fasta)

            # Cnfigure result TSV lines
            result_tsv_lines = parse_align_results_xml(align_xml_text,
                packet["qual"], taxonomy_path)

            # Write the result to tsv
            write_classification(result_tsv_lines, tsv_res_path)
        # end for

        with print_lock:
            printl(logfile_path, "\r{} - File '{}' is processed.".format(getwt(), os.path.basename(fq_fa_path)))
            printn("Working...")
        # end with
    # end for
    remove_tmp_files( os.path.join(queries_tmp_dir, "query{}_tmp.fasta".format(os.getpid())) )
# end def process_paral


def process(fq_fa_list, n_thr, packet_size, tax_annot_res_dir,
    blast_algorithm, use_index, logfile_path):
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
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    print_lock = mp.Lock() # lock for printing

    pool = mp.Pool(n_thr, initializer=init_process,
        initargs=(print_lock, packet_size,
        tax_annot_res_dir, blast_algorithm,
        use_index, logfile_path))

    pool.starmap(process_paral, [ (fq_fa_sublist,) for fq_fa_sublist in spread_files_equally(fq_fa_list, n_thr) ])

    # Reaping zombies
    pool.close()
    pool.join()
# end def process
