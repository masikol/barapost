# -*- coding: utf-8 -*-
# Module defines functions necessary for sorting FASTA and FASTQ files in single thread.

import os
import sys

from src.sorter_modules.sorter_spec import *

from src.fmt_readID import fmt_read_id
from src.platform import platf_depend_exit
from src.printlog import printl, printn, getwt, err_fmt
from src.filesystem import get_curr_res_dpath, is_fastq

from src.sorter_modules.fastq_records import fastq_records
from src.sorter_modules.fasta_records import fasta_records


def write_fastq_record(sorted_file, fastq_record):
    """
    :param sorted_path: path to file, in which data from fastq_record is oing to be written;
    :type sorted_path: str;
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of FASTQ;
    :type fastq_record: dict<str: str>;
    """

    sorted_file.write(fastq_record["seq_id"]+'\n')
    sorted_file.write(fastq_record["seq"]+'\n')
    sorted_file.write(fastq_record["opt_id"]+'\n')
    sorted_file.write(fastq_record["qual_line"]+'\n')
# end def write_fastq_record


def write_fasta_record(sorted_file, fasta_record):
    """
    :param sorted_file: file, which data from fasta_record is written in
    :type sorted_file: _io.TextIOWrapper
    :param fasta_record: dict of 2 elements. Elements are four corresponding lines of FASTA
    :type fasta_record: dict<str: str>
    """
    sorted_file.write(fasta_record["seq_id"]+'\n')
    sorted_file.write(fasta_record["seq"]+'\n')

# end def write_fasta_record


def sort_fastqa_file(fq_fa_path, tax_annot_res_dir, sens,
        min_qual, min_qlen, logfile_path):
    """
    Function for single-thread sorting FASTQ and FASTA files.

    :param fq_fa_path: path to FASTQ (of FASTA) file meant to be processed;
    :type fq_fa_path: str;
    :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    :type tax_annot_res_dir: str;
    :param sens: sorting sensitivity;
    :type sens: str;
    :param min_qual: minimum quality to keep;
    :type min_qual: float;
    :param min_qlen: minimmum sequence length to keep;
    :type min_qlen: int (or None, if this feature is disabled);
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    outdir_path = os.path.dirname(logfile_path)
    minlen_fmt_str = "_len_less_{}".format(min_qlen) if not min_qlen is None else ""

    seqs_pass = 0
    seqs_fail = 0
    srt_file_dict = dict() # dict containing file objects of existing output files

    new_dpath = get_curr_res_dpath(fq_fa_path, tax_annot_res_dir)
    tsv_res_fpath = get_res_tsv_fpath(new_dpath)
    resfile_lines = configure_resfile_lines(tsv_res_fpath, sens)

    # Configure generator, write function and path to trash file
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

    for fastq_rec in seq_records_generator(fq_fa_path):

        read_name = sys.intern(fmt_read_id(fastq_rec["seq_id"])[1:]) # get ID of the sequence

        try:
            hit_names, quality, q_len = resfile_lines[read_name]  # find hit corresponding to this sequence
        except KeyError:
            printl(logfile_path, err_fmt("""read '{}' not found in TSV file containing taxonomic annotation.
This TSV file: '{}'""".format(read_name, tsv_res_fpath)))
            printl(logfile_path, "Make sure that this read has been already processed by 'prober.py' and 'barapost.py'.")
            platf_depend_exit(1)
        # end try

        # If read is found in TSV file:
        q_len = SeqLength(q_len)
        if q_len < min_qlen or (quality != '-' and quality < min_qual):
            # Place this sequence to trash file
            if trash_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, trash_fpath)
            # end if
            write_fun(srt_file_dict[trash_fpath], fastq_rec) # write current read to sorted file
            seqs_fail += 1
        else:
            for hit_name in hit_names.split("&&"): # there can be multiple hits for single query sequence
                # Get name of result FASTQ file to write this read in
                sorted_file_path = os.path.join(outdir_path, "{}.fast{}".format(hit_name,
                    'q' if is_fastq(fq_fa_path) else 'a'))
                if sorted_file_path not in srt_file_dict.keys():
                    srt_file_dict = update_file_dict(srt_file_dict, sorted_file_path)
                # end if
                write_fun(srt_file_dict[sorted_file_path], fastq_rec) # write current read to sorted file
            # end for
            seqs_pass += 1
        # end if
    # end for

    # Close all sorted files
    for file_obj in srt_file_dict.values():
        file_obj.close()
    # end for

    printl(logfile_path, "\r{} - File '{}' is sorted.".format(getwt(), os.path.basename(fq_fa_path)))
    printn(" Working...")

    return (seqs_pass, seqs_fail)
# end def sort_fastqa_file
