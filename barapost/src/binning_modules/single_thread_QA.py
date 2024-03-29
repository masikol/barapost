# -*- coding: utf-8 -*-
# Module defines functions necessary for binning FASTA and FASTQ files in single thread.

import os
import sys
import logging

from src.binning_modules.binning_spec import get_res_tsv_fpath, configure_resfile_lines

from src.fmt_read_id import fmt_read_id
from src.platform import platf_depend_exit
from src.printlog import printlog_error, printlog_error_time
from src.filesystem import get_curr_res_dpath, is_fastq

from src.binning_modules.fastq_records import fastq_records
from src.binning_modules.fasta_records import fasta_records

from src.binning_modules.filters import get_QL_filter, get_QL_trash_fpath
from src.binning_modules.filters import get_align_filter, get_align_trash_fpath
from src.binning_modules.filters import get_classif_not_found_fpath


def write_fastq_record(binned_file, fastq_record):
    # :param binned_file: file instance, in which data from fastq_record is oing to be written;
    # :type binned_file: _oi.TextIOWrapper;
    # :param fastq_record: dict of 4 elements. Elements are four corresponding lines of FASTQ;
    # :type fastq_record: dict<str: str>;

    if not binned_file is None:
        binned_file.write(fastq_record["seq_id"]+'\n')
        binned_file.write(fastq_record["seq"]+'\n')
        binned_file.write(fastq_record["opt_id"]+'\n')
        binned_file.write(fastq_record["qual_line"]+'\n')
    # end if
# end def write_fastq_record


def write_fasta_record(binned_file, fasta_record):
    # :param binned_file: file, which data from fasta_record is written in
    # :type binned_file: _io.TextIOWrapper
    # :param fasta_record: dict of 2 elements. Elements are four corresponding lines of FASTA
    # :type fasta_record: dict<str: str>

    if not binned_file is None:
        binned_file.write(fasta_record["seq_id"]+'\n')
        binned_file.write(fasta_record["seq"]+'\n')
    # end if
# end def write_fasta_record


def update_file_dict(srt_file_dict, new_fpath):
    try:
        if not new_fpath is None:
            srt_file_dict[sys.intern(new_fpath)] = open(new_fpath, 'a')
        else:
            srt_file_dict[new_fpath] = None # handle no_trash
        # end if
    except OSError as oserr:
        printlog_error_time("Error occured while opening one of result files")
        printlog_error("Errorneous file: `{}`".format(new_fpath))
        printlog_error( str(oserr) )
        platf_depend_exit(1)
    # end try
    return srt_file_dict
# end def update_file_dict


def bin_fastqa_file(fq_fa_path, tax_annot_res_dir, sens,
        min_qual, min_qlen, min_pident, min_coverage, no_trash):
    # Function for single-thread binning FASTQ and FASTA files.
    #
    # :param fq_fa_path: path to FASTQ (of FASTA) file meant to be processed;
    # :type fq_fa_path: str;
    # :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    # :type tax_annot_res_dir: str;
    # :param sens: binning sensitivity;
    # :type sens: str;
    # :param min_qual: threshold for quality filter;
    # :type min_qual: float;
    # :param min_qlen: threshold for length filter;
    # :type min_qlen: int (or None, if this filter is disabled);
    # :param min_pident: threshold for alignment identity filter;
    # :type min_pident: float (or None, if this filter is disabled);
    # :param min_coverage: threshold for alignment coverage filter;
    # :type min_coverage: float (or None, if this filter is disabled);
    # :param no_trash: loical value. True if user does NOT want to output trash files;
    # :type no_trash: bool;

    outdir_path = os.path.dirname(logging.getLoggerClass().root.handlers[0].baseFilename)

    seqs_pass = 0 # counter for sequences, which pass filters
    QL_seqs_fail = 0 # counter for too short or too low-quality sequences
    align_seqs_fail = 0 # counter for sequences, which align to their best hit with too low identity or coverage

    srt_file_dict = dict() # dict containing file objects of existing output files

    new_dpath = get_curr_res_dpath(fq_fa_path, tax_annot_res_dir)
    tsv_res_fpath = get_res_tsv_fpath(new_dpath)
    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy", "taxonomy.tsv")
    resfile_lines = configure_resfile_lines(tsv_res_fpath, sens, taxonomy_path)

    # Configure generator, write function and path to trash file
    if is_fastq(fq_fa_path):
        seq_records_generator = fastq_records
        write_fun =  write_fastq_record
    else:
        seq_records_generator = fasta_records
        write_fun = write_fasta_record
    # end if

    # Configure path to "classification not found" file
    classif_not_found_fpath = get_classif_not_found_fpath(fq_fa_path, outdir_path)

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

    for fastq_rec in seq_records_generator(fq_fa_path):

        read_name = sys.intern(fmt_read_id(fastq_rec["seq_id"])[1:]) # get ID of the sequence

        try:
            hit_names, *vals_to_filter = resfile_lines[read_name]  # find hit corresponding to this sequence
        except KeyError:
            # Place this sequence into the "classification not found" file
            if classif_not_found_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, classif_not_found_fpath)
            # end if
            write_fun(srt_file_dict[classif_not_found_fpath], fastq_rec) # write current read to binned file
            continue
        # end try

        # Apply filters
        if not QL_filter(vals_to_filter):
            QL_seqs_fail += 1
            # Place this sequence to QL trash file
            if QL_trash_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, QL_trash_fpath)
            # end if
            write_fun(srt_file_dict[QL_trash_fpath], fastq_rec) # write current read to binned file

        elif not align_filter(vals_to_filter):
            align_seqs_fail += 1
            # Place this sequence to align_trash file
            if align_trash_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, align_trash_fpath)
            # end if
            write_fun(srt_file_dict[align_trash_fpath], fastq_rec) # write current read to binned file

        else:
            for hit_name in hit_names.split("&&"): # there can be multiple hits for single query sequence
                # Get name of result FASTQ file to write this read in
                binned_file_path = os.path.join(outdir_path, "{}.fast{}".format(hit_name,
                    'q' if is_fastq(fq_fa_path) else 'a'))
                if binned_file_path not in srt_file_dict.keys():
                    srt_file_dict = update_file_dict(srt_file_dict, binned_file_path)
                # end if
                write_fun(srt_file_dict[binned_file_path], fastq_rec) # write current read to binned file
            # end for
            seqs_pass += 1
        # end if
    # end for

    # Close all binned files
    for file_obj in filter(lambda x: not x is None, srt_file_dict.values()):
        file_obj.close()
    # end for

    return (seqs_pass, QL_seqs_fail, align_seqs_fail)
# end def bin_fastqa_file
