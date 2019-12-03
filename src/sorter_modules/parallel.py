# -*- coding: utf-8 -*-

# If we process files in parallel, we need to open output files each time we need to write to it,
#    because different processed cannot correctly write to the same file and have their own file descriptors:
#    it will produce broken files.
# Thus we will open these files each time and descriptors will be up-to-date.


from gzip import open as open_as_gzip
from re import search as re_search
from glob import glob
import os


is_gzipped = lambda f: True if f.endswith(".gz") else False
OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line,   # format text line
    lambda line: line.decode("utf-8")  # format gzipped line
)


def printl(text=""):
    """
    Function for printing text to console and to log file.
    """
    print(text)
    logfile.write(str(text).strip('\r') + '\n')
    logfile.flush()
# end def printl


from src.sorter_modules.common import *


def write_fastq_record(sorted_path, fastq_record):
    """
    :param sorted_path: path to file, in which data from fastq_record is oing to be written;
    :type sorted_path: str;
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of FASTQ;
    :type fastq_record: dict<str: str>;
    """
    try:
        with open(sorted_path, 'a') as sorted_file:

            sorted_file.write(fastq_record["seq_id"])
            sorted_file.write(fastq_record["seq"])
            sorted_file.write(fastq_record["opt_id"])
            sorted_file.write(fastq_record["qual_line"])
        # end with
    except OSError as err:
        printl(err_fmt( str(err) ))
        printl("File: '{}'".format(sorted_path))
        platf_depend_exit(0)
    # end try
# end def write_fastq_record


def write_fasta_record(sorted_path, fasta_record):
    """
    :param sorted_path: path to file, in which data from fastq_record is oing to be written;
    :type sorted_path: str;
    :param fasta_record: dict of 2 elements. Elements are four corresponding lines of FASTA;
    :type fasta_record: dict<str: str>;
    """
    try:
        with open(sorted_path, 'a') as sorted_file:
            sorted_file.write(fasta_record["seq_id"])
            sorted_file.write(fasta_record["seq"])
        # end with
    except OSError as err:
        printl(err_fmt( str(err) ))
        printl("File: '{}'".format(sorted_path))
        platf_depend_exit(0)
    # end try
# end def write_fasta_record


def spread_files_equally(fq_fa_list, n_thr):
    """
    Function distributes files among processes equally.

    :param fq_fa_list: list of paths to files meant to be processed:
    :type fq_fa_list: list<str>;
    :param n_thr: number of therads to launch;
    :type n_thr: int;
    """

    sublist_size = len(fq_fa_list) // n_thr

    # Processes [0, (n_thr-1)] will obtain equally 'sublist_size' files:
    start_pos = 0
    for i in range(n_thr - 1):
        yield fq_fa_list[start_pos : start_pos+sublist_size]
        start_pos += sublist_size
    # end for

    # Give the rest of data to the last unlucky process:
    yield fq_fa_list[start_pos :]
# end def spread_files_equally


def init_paral_sorting(write_lock_buff, inc_val_buff, inc_val_lock_buff):

    global write_lock
    write_lock = write_lock_buff

    global inc_val
    inc_val = inc_val_buff

    global inc_val_lock
    inc_val_lock = inc_val_lock_buff
# end def init_paral_sorting


def sort_fastqa_file(fq_fa_path):
    """
    Function for parallel sorting FASTQ and FASTA files.

    :param fq_fa_path: path to FASTQ (of FASTA) file meant to be processed;
    :type fq_fa_path: str;
    """

    seqs_pass = 0
    seqs_fail = 0

    new_dpath = get_curr_res_dpath(fq_fa_path, tax_annot_res_dir)
    tsv_res_fpath = get_res_tsv_fpath(new_dpath)
    resfile_lines = configure_resfile_lines(tsv_res_fpath, sens)

    # Configure path to trash file
    if is_fastq(fq_fa_path):
        seq_records_generator = fastq_records
        write_fun =  write_fastq_record
        trash_fpath = os.path.join(outdir_path, "qual_less_Q{}{}.fastq".format(int(min_ph33_qual),
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
                hit_name, ph33_qual, q_len = resfile_lines[read_name] # find hit corresponding to this sequence
            except KeyError:
                printl(err_fmt("""read '{}' not found in TSV file containing taxonomic annotation.
This TSV file: '{}'""".format(read_name, tsv_res_fpath)))
                printl("Make sure that this read has been already processed by 'prober.py' and 'barapost.py'.")
                platf_depend_exit(1)
            # end try

            # If read is found in TSV file:
            q_len = SeqLength(q_len)
            if q_len < min_qlen or (ph33_qual != '-' and ph33_qual < min_ph33_qual):
                # Place this sequence to trash file
                to_write[read_name] = (fastqa_rec, trash_fpath)
                seqs_fail += 1
            else:
                # Get name of result FASTQ file to write this read in
                sorted_file_path = os.path.join(outdir_path, "{}.fast{}".format(hit_name,
                    'q' if is_fastq(fq_fa_path) else 'a'))
                to_write[read_name] = (fastqa_rec, sorted_file_path)
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

    return (seqs_pass, seqs_fail)
# end def sort_fastqa_file