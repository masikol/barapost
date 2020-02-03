# -*- coding: utf-8 -*-

import os
import multiprocessing as mp
from re import search as re_search

from src.printlog import getwt, printl, printn
from src.filesystem import get_curr_res_dpath, create_result_directory, remove_tmp_files, is_fastq
from src.write_classification import write_classification
from src.filesystem import OPEN_FUNCS, is_gzipped

from src.barapost_modules.barapost_spec import look_around, launch_blastn, parse_align_results_xml

from src.barapost_modules.fasta_packets_from_str import fasta_packets_from_str
from src.fasta import fasta_packets
from src.fastq import fastq_packets


def init_proc_single_file_in_paral(print_lock_buff, write_lock_buff, packet_size_buff,
    tax_annot_res_dir_buff, blast_algorithm_buff, use_index_buff, logfile_path_buff):
    """
    Function that initializes global variables that all processes shoud have access to.
    This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.
    Function works when 'few_files'-parallel mode is running.

    :param print_lock_buff: lock that synchronizes printing to the console;
    :type print_lock_buff: multiprocessing.Lock;
    :param write_lock_buff: lock that synchronizes wriiting to result file;
    :type write_lock_buff: multiprocessing.Lock;
    """
    global print_lock
    print_lock = print_lock_buff

    global write_lock
    write_lock = write_lock_buff

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
# end def init_proc_single_file_in_paral


def process_part_of_file(data, tsv_res_path):
    """
    Function preforms processing part of file in 'few_files'-parallel mode.

    :param data: FASTA data meant to be processed;
    :type data: str;
    :param tsv_res_path: path to result TSV file;
    :type tsv_res_path: str;
    """

    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy","taxonomy")
    queries_tmp_dir = os.path.join(tax_annot_res_dir, "queries-tmp")
    local_fasta = os.path.join(tax_annot_res_dir, "local_database", "local_seq_set.fasta")

    for packet in fasta_packets_from_str(data, packet_size, num_done_seqs=0):

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
        with write_lock:
            write_classification(result_tsv_lines, tsv_res_path)
        # end with
    # end for

    remove_tmp_files( os.path.join(queries_tmp_dir, "query{}_tmp.fasta".format(os.getpid())) )
# end def process_part_of_file


def process(fq_fa_path, n_thr, packet_size, tax_annot_res_dir,
            blast_algorithm, use_index, logfile_path):
    """
    Function preforms "few_files"-parallel mode.

    :param fq_fa_path: path to FASTA or FASTQ file meant to be processed;
    :type fq_fa_path: str;
    :param i: number of this file;
    :type i: int;
    """

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
        tsv_res_path = "{}.tsv".format(os.path.join(new_dpath,
            "classification")) # form result tsv file path
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
        printl(logfile_path, "\rFile '{}' have been already completely processed.".format(fq_fa_path))
        printl(logfile_path, "Omitting it.")
        printn("  Working...")
        return
    # end if

    file_part_size = num_seqs // n_thr
    if num_seqs // n_thr != 0:
        file_part_size += 1
    # end if

    print_lock = mp.Lock() # lock for printing to the console
    write_lock = mp.Lock() # lock for writing to the result file

    pool = mp.Pool(n_thr, initializer=init_proc_single_file_in_paral,
        initargs=(print_lock, write_lock, packet_size, tax_annot_res_dir,
            blast_algorithm, use_index, logfile_path))

    pool.starmap(process_part_of_file, [(file_part["fasta"],
        tsv_res_path) for file_part in packet_generator(fq_fa_path,
        file_part_size, # part of file instead of actually packet size (piece of data to Process)
        num_done_seqs)])

    # Reaping zombies
    pool.close()
    pool.join()

    printl(logfile_path, "\r{} - File '{}' is processed.".format(getwt(), os.path.basename(fq_fa_path)))
    printn("  Working...")
# end def process
