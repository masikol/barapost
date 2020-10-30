# -*- coding: utf-8 -*-
# Module defines functions necessary for binning FASTA and FASTQ files in parallel.

# If we process files in parallel, we need to open output files each time we need to write to it,
#    because different processed cannot correctly write to the same file and have their own file descriptors:
#    it will produce broken files.
# Thus we will open these files each time and descriptors will be up-to-date.

import os
import sys
import multiprocessing as mp

from src.binning_modules.binning_spec import *

from src.fmt_readID import fmt_read_id
from src.platform import platf_depend_exit
from src.printlog import printl, printn, getwt, err_fmt
from src.filesystem import get_curr_res_dpath, is_fastq

from src.binning_modules.fastq_records import fastq_records
from src.binning_modules.fasta_records import fasta_records

from src.binning_modules.filters import get_QL_filter, get_QL_trash_fpath
from src.binning_modules.filters import get_align_filter, get_align_trash_fpath


def write_fastq_record(binned_fpath, fastq_record):
    """
    :param binned_fpath: path to file, in which data from fastq_record is oing to be written;
    :type binned_fpath: str;
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of FASTQ;
    :type fastq_record: dict<str: str>;
    """

    if not binned_fpath is None:
        with open(binned_fpath, 'a') as binned_file:
            binned_file.write(fastq_record["seq_id"]+'\n')
            binned_file.write(fastq_record["seq"]+'\n')
            binned_file.write(fastq_record["opt_id"]+'\n')
            binned_file.write(fastq_record["qual_line"]+'\n')
        # end with
    # end if
# end def write_fastq_record


def write_fasta_record(binned_fpath, fasta_record):
    """
    :param binned_file: file, which data from fasta_record is written in
    :type binned_file: _io.TextIOWrapper
    :param fasta_record: dict of 2 elements. Elements are four corresponding lines of FASTA
    :type fasta_record: dict<str: str>
    """

    if not binned_fpath is None:
        with open(binned_fpath, 'a') as binned_file:
            binned_file.write(fasta_record["seq_id"]+'\n')
            binned_file.write(fasta_record["seq"]+'\n')
        # end with
    # end if
# end def write_fasta_record


def init_paral_binning(print_lock_buff, write_lock_buff):
    """
    Function initializes global locks for parallel binning of fasta and fastq files.

    :param print_lock_buff: lock for printing to console;
    :type print_lock_buff: multiprocessing.Lock;
    :param write_lock_buff: lock for printing writing to binned files;
    :type write_lock_buff: multiprocessing.Lock;
    """

    global print_lock
    print_lock = print_lock_buff
    
    global write_lock
    write_lock = write_lock_buff

# end def init_paral_binning


def bin_fastqa_file(fq_fa_lst, tax_annot_res_dir, sens, n_thr,
    min_qual, min_qlen, min_pident, min_coverage, no_trash, logfile_path):
    """
    Function for parallel binning FASTQ and FASTA files.
    Actually bins multiple files.

    :param fq_fa_lst: lsit of paths to FASTQ (of FASTA) file meant to be processed;
    :type fq_fa_lst: list<str>;
    :param min_qual: threshold for quality filter;
    :type min_qual: float;
    :param min_qlen: threshold for length filter;
    :type min_qlen: int (or None, if this filter is disabled);
    :param min_pident: threshold for alignment identity filter;
    :type min_pident: float (or None, if this filter is disabled);
    :param min_coverage: threshold for alignment coverage filter;
    :type min_coverage: float (or None, if this filter is disabled);
    :param no_trash: loical value. True if user does NOT want to output trash files;
    :type no_trash: bool;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    outdir_path = os.path.dirname(logfile_path)

    seqs_pass = 0 # counter for sequences, which pass filters
    QL_seqs_fail = 0 # counter for too short or too low-quality sequences
    align_seqs_fail = 0 # counter for sequences, which align to their best hit with too low identity or coverage

    for fq_fa_path in fq_fa_lst:

        new_dpath = get_curr_res_dpath(fq_fa_path, tax_annot_res_dir)
        tsv_res_fpath = get_res_tsv_fpath(new_dpath)
        resfile_lines = configure_resfile_lines(tsv_res_fpath, sens,
            os.path.join(tax_annot_res_dir, "taxonomy", "taxonomy"), logfile_path)

        # Configure path to trash file
        if is_fastq(fq_fa_path):
            seq_records_generator = fastq_records
            write_fun =  write_fastq_record
        else:
            seq_records_generator = fasta_records
            write_fun = write_fasta_record
        # end if

        # Make filter for quality and length
        QL_filter = get_QL_filter(fq_fa_path, min_qual, min_qlen)
        # Configure path to trash file
        if not no_trash:
            QL_trash_fpath = get_QL_trash_fpath(fq_fa_path, outdir_path, min_qual, min_qlen,)
        else:
            QL_trash_fpath = None
        # end if

        # Make filter for identity and coverage
        align_filter = get_align_filter(min_pident, min_coverage)
        # Configure path to this trash file
        if not no_trash:
            align_trash_fpath = get_align_trash_fpath(fq_fa_path, outdir_path, min_pident, min_coverage)
        else:
            align_trash_fpath = None
        # end if

        # Create an iterator that will yield records
        seq_records_iterator = iter(seq_records_generator(fq_fa_path))
        # Dict for storing batches of sequences meant to be written to output files:
        to_write = dict()
        stop = False # for outer while-loop

        while not stop:

            # Extract batch of records of 'n_thr' size and find their destination paths:
            for _ in range(n_thr):

                try:
                    fastqa_rec = next(seq_records_iterator)
                except StopIteration:
                    stop = True # for outer while-loop
                    break
                # end try

                read_name = sys.intern(fmt_read_id(fastqa_rec["seq_id"])[1:]) # get ID of the sequence

                try:
                    hit_names, *vals_to_filter = resfile_lines[read_name]  # find hit corresponding to this sequence
                except KeyError:
                    printl(logfile_path, err_fmt("""read '{}' not found in TSV file containing taxonomic annotation.
    This TSV file: '{}'""".format(read_name, tsv_res_fpath)))
                    printl(logfile_path, "Make sure that this read has been already processed by 'barapost-prober.py' and 'barapost-local.py'.")
                    platf_depend_exit(1)
                # end try

                # If read is found in TSV file:
                if not QL_filter(vals_to_filter):
                    # Place this sequence to QL trash file
                    to_write[read_name] = (fastqa_rec, QL_trash_fpath)
                    QL_seqs_fail += 1
                elif not align_filter(vals_to_filter):
                    # Place this sequence to QL trash file
                    to_write[read_name] = (fastqa_rec, align_trash_fpath)
                    align_seqs_fail += 1
                else:
                    for hit_name in hit_names.split("&&"):
                        # Get name of result FASTQ file to write this read in
                        binned_file_path = os.path.join(outdir_path, "{}.fast{}".format(hit_name,
                            'q' if is_fastq(fq_fa_path) else 'a'))
                        to_write[read_name] = (fastqa_rec, binned_file_path)
                    # end for
                    seqs_pass += 1
                # end if
            # end for

            # Write batch of records to output files:
            with write_lock:
                for record, fpath in to_write.values():
                    write_fun(fpath, record)
                # end for
            # end with
            to_write.clear()
        # end while

        with write_lock:
            # Write the rest of 'uneven' data to output files:
            if len(to_write) != 0:
                for record, fpath in to_write.values():
                    write_fun(fpath, record)
                # end for
            # end if
            printl(logfile_path, "\r{} - File '{}' is binned.".format(getwt(), os.path.basename(fq_fa_path)))
            printn(" Working...")
        # end with
    # end for

    return (seqs_pass, QL_seqs_fail, align_seqs_fail)
# end def bin_fastqa_file
