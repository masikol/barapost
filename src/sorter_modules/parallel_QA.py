# -*- coding: utf-8 -*-

# If we process files in parallel, we need to open output files each time we need to write to it,
#    because different processed cannot correctly write to the same file and have their own file descriptors:
#    it will produce broken files.
# Thus we will open these files each time and descriptors will be up-to-date.

from src.sorter_modules.sorter_spec import *

from src.printlog import printl, printn, getwt, err_fmt
from src.platform import platf_depend_exit
from src.filesystem import get_curr_res_dpath, is_fastq
from src.fmt_readID import fmt_read_id

from src.sorter_modules.fastq_records import fastq_records
from src.sorter_modules.fasta_records import fasta_records

import os


def write_fastq_record(sorted_path, fastq_record):
    """
    :param sorted_path: path to file, in which data from fastq_record is oing to be written;
    :type sorted_path: str;
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of FASTQ;
    :type fastq_record: dict<str: str>;
    """

    with open(sorted_path, 'a') as sorted_file:

        sorted_file.write(fastq_record["seq_id"]+'\n')
        sorted_file.write(fastq_record["seq"]+'\n')
        sorted_file.write(fastq_record["opt_id"]+'\n')
        sorted_file.write(fastq_record["qual_line"]+'\n')
    # end with
# end def write_fastq_record


def write_fasta_record(sorted_fpath, fasta_record):
    """
    :param sorted_file: file, which data from fasta_record is written in
    :type sorted_file: _io.TextIOWrapper
    :param fasta_record: dict of 2 elements. Elements are four corresponding lines of FASTA
    :type fasta_record: dict<str: str>
    """

    with open(sorted_fpath, 'a') as sorted_file:
        sorted_file.write(fasta_record["seq_id"]+'\n')
        sorted_file.write(fasta_record["seq"]+'\n')
    # end with
# end def write_fasta_record


def init_paral_sorting(print_lock_buff, write_lock_buff):

    global print_lock
    print_lock = print_lock_buff
    
    global write_lock
    write_lock = write_lock_buff

# end def init_paral_sorting


def sort_fastqa_file(fq_fa_lst, tax_annot_res_dir, sens, n_thr,
    min_qual, min_qlen, logfile_path):
    """
    Function for parallel sorting FASTQ and FASTA files.
    Actually sorts multiple files.

    :param fq_fa_lst: lsit of paths to FASTQ (of FASTA) file meant to be processed;
    :type fq_fa_lst: list<str>;
    """

    outdir_path = os.path.dirname(logfile_path)
    minlen_fmt_str = "_len_less_{}".format(min_qlen) if not min_qlen is None else ""

    seqs_pass = 0
    seqs_fail = 0

    for fq_fa_path in fq_fa_lst:

        new_dpath = get_curr_res_dpath(fq_fa_path, tax_annot_res_dir)
        tsv_res_fpath = get_res_tsv_fpath(new_dpath)
        resfile_lines = configure_resfile_lines(tsv_res_fpath, sens)

        # Configure path to trash file
        if is_fastq(fq_fa_path):
            seq_records_generator = fastq_records
            write_fun =  write_fastq_record
            trash_fpath = os.path.join(outdir_path, "qual_less_Q{}{}.fastq".format(int(min_qual),
                minlen_fmt_str))
        else:
            seq_records_generator = fasta_records
            write_fun = write_fasta_record
            trash_fpath = os.path.join(outdir_path, "len_less_{}.fasta".format(min_qlen))
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

                read_name = sys.intern(fmt_read_id(fastqa_rec["seq_id"])) # get ID of the sequence

                try:
                    hit_names, ph33_qual, q_len = resfile_lines[read_name[1:]] # find hit corresponding to this sequence
                except KeyError:
                    printl(logfile_path, err_fmt("""read '{}' not found in TSV file containing taxonomic annotation.
    This TSV file: '{}'""".format(read_name, tsv_res_fpath)))
                    printl(logfile_path, "Make sure that this read has been already processed by 'prober.py' and 'barapost.py'.")
                    platf_depend_exit(1)
                # end try

                # If read is found in TSV file:
                q_len = SeqLength(q_len)
                if q_len < min_qlen or (ph33_qual != '-' and ph33_qual < min_qual):
                    # Place this sequence to trash file
                    to_write[read_name] = (fastqa_rec, trash_fpath)
                    seqs_fail += 1
                else:
                    for hit_name in hit_names.split("&&"):
                        # Get name of result FASTQ file to write this read in
                        sorted_file_path = os.path.join(outdir_path, "{}.fast{}".format(hit_name,
                            'q' if is_fastq(fq_fa_path) else 'a'))
                        to_write[read_name] = (fastqa_rec, sorted_file_path)
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

        # Write the rest of 'uneven' data to output files:
        if len(to_write) != 0:
            with write_lock:
                for record, fpath in to_write.values():
                    write_fun(fpath, record)
                # end for
            # end with
        # end if

        printl(logfile_path, "\r{} - File '{}' is sorted.".format(getwt(), os.path.basename(fq_fa_path)))
        printn(" Working...")
    # end for


    return (seqs_pass, seqs_fail)
# end def sort_fastqa_file
