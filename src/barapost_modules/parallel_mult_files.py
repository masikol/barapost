# -*- coding: utf-8  -*-

from src.printlog import getwt, printl, printn
from src.spread_files_equally import spread_files_equally
from src.filesystem import get_curr_res_dpath, create_result_directory, remove_tmp_files, is_fastq
from src.write_classification import write_classification

from src.fasta import fasta_packets
from src.fastq import fastq_packets

from src.barapost_modules import look_around, launch_blastn, parse_align_results_xml

import os
import multiprocessing as mp
from re import search as re_search

def process(fq_fa_list, n_thr, packet_size, tax_annot_res_dir,
    blast_algorithm, use_index, logfile_path):

    print_lock = mp.Lock() # lock for printing

    pool = mp.Pool(n_thr, initializer=init_process, initargs=(print_lock, packet_size, tax_annot_res_dir,
        blast_algorithm, use_index, logfile_path))
    pool.starmap(process, [ (fq_fa_sublist,) for fq_fa_sublist in spread_files_equally(fq_fa_list, n_thr) ])

    # Reaping zombies
    pool.close()
    pool.join()
# end def process


def init_process(print_lock_buff, packet_size_buff, tax_annot_res_dir_buff, blast_algorithm_buff, use_index_buff):
    """
    Function that initializes global variables that all processes shoud have access to.
    This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.
    Function works when 'many_files'-parallel mode is running.

    :param print_lock_buff: lock that synchronizes printing to the console;
    :type print_lock_buff: multiprocessing.Lock;
    """

    global print_lock
    print_lock = print_lock_buff

    global tax_annot_res_dir
    tax_annot_res_dir = tax_annot_res_dir_buff

    global packet_size
    packet_size = packet_size_buff

    global blast_algorithm
    blast_algorithm = blast_algorithm_buff

    global use_index
    use_index = use_index_buff

# end def init_proc_many_files


def process_paral(fq_fa_list):
    """
    Function performs 'many_files'-parallel mode of single-thread mode.
    They differ only in ptinting to the console.

    :param fq_fa_list: list of paths to FASTA and FASTQ files meant to be processed;
    :type fq_fa_list: list<str>;
    :param parallel: flag indicating if parallel mode if enabled.
        Influences only on printing to the console;
    :type parallel: bool;
    """

    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy","taxonomy")
    queries_tmp_dir = os.path.join(tax_annot_res_dir, "queries-tmp")

    # Iterate over source FASTQ and FASTA files
    for i, fq_fa_path in enumerate(fq_fa_list):

        # Create the result directory with the name of FASTQ of FASTA file being processed:
        new_dpath = create_result_directory(fq_fa_path, tax_annot_res_dir)

        # "hname" means human readable name (i.e. without file path and extention)
        infile_hname = os.path.basename(fq_fa_path)
        infile_hname = re_search(r"(.+)\.(m)?f(ast)?(a|q)(\.gz)?$", infile_hname).group(1)

        # Look around and ckeck if there are results of previous runs of this script
        # If 'look_around' is None -- there is no data from previous run
        previous_data = look_around(new_dpath, fq_fa_path, blast_algorithm, logfile_path)

        if previous_data is None: # If there is no data from previous run
            num_done_reads = 0 # number of successfully processed sequences
            tsv_res_path = "{}.tsv".format(os.path.join(new_dpath,
                "classification")) # form result tsv file path
        else: # if there is data from previous run
            num_done_reads = previous_data["n_done_reads"] # get number of successfully processed sequences
            tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        # end if

        packet_generator = fastq_packets if is_fastq(fq_fa_path) else fasta_packets

        for packet in packet_generator(fq_fa_path, packet_size, num_done_reads):

            if packet["fasta"] == "":
                with print_lock:
                    printl(logfile_path, "\nFile '{}' has been already completely processed.".format(fq_fa_path))
                    printl(logfile_path, "Omitting it.")
                # end with
                continue
            # end if

            # Align the packet
            align_xml_text = launch_blastn(packet["fasta"], blast_algorithm,
                use_index, queries_tmp_dir, local_fasta)

            # Get result tsv lines
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
