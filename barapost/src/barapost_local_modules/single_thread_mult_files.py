# -*- coding: utf-8  -*-
# This module defines functions necessary for barapost-local.py to perform single-thread work.

import os
import re
import sys

from src.printlog import printlog_info, printlog_info_time, printn, printlog_warning
from src.write_classification import write_classification
from src.filesystem import create_result_directory
from src.filesystem import remove_tmp_files, is_fastq, is_gzipped, OPEN_FUNCS, FORMATTING_FUNCS

from src.fasta import fasta_packets
from src.fastq import fastq_packets

from src.barapost_local_modules.barapost_spec import look_around, launch_blastn, parse_align_results_xml


def process(fq_fa_list, packet_size, tax_annot_res_dir, blast_algorithm, use_index, db_path):
    # Function launches parallel processing in "many-files" mode by barapost-local.py.
    #
    # :param fq_fa_list: list of paths to files meant to be processed;
    # :type fq_fa_list: list<str>;
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

    nfiles = len(fq_fa_list)

    # Iterate over source FASTQ and FASTA files
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
            tsv_res_path = "{}.tsv".format(os.path.join(new_dpath,
                "classification")) # form result tsv file path
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

        if num_seqs == num_done_seqs:
            sys.stdout.write('\r')
            printlog_info_time("File #{}/{} (`{}`) has been already completely processed."\
                .format(i+1, nfiles, fq_fa_path))
            printlog_info("Omitting it.")
            printn("Working...")
            continue
        # end if

        for packet in packet_generator(fq_fa_path, packet_size, num_done_seqs):

            # Align the packet
            align_xml_text = launch_blastn(packet["fasta"], blast_algorithm,
                use_index, queries_tmp_dir, db_path)

            # Get result tsv lines
            result_tsv_lines = parse_align_results_xml(align_xml_text,
                packet["qual"])

            # Write the result to tsv
            write_classification(result_tsv_lines, tsv_res_path)
        # end for

        sys.stdout.write('\r')
        printlog_info_time("File #{}/{} (`{}`) is processed."\
            .format(i+1, nfiles, os.path.basename(fq_fa_path)))
        printn("Working...")
    # end for

    query_fpath = os.path.join(queries_tmp_dir, "query{}_tmp.fasta".format(os.getpid()))
    remove_tmp_files(query_fpath)
# end def process
