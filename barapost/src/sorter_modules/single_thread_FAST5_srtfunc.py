# -*- coding: utf-8 -*-
# Module defines functions necessary for sorting FAST5 files "directly": without index.

import h5py

import os
import sys
from glob import glob

from src.sorter_modules.sorter_spec import *
from src.sorter_modules.fast5 import update_file_dict

from src.platform import platf_depend_exit
from src.printlog import printl, printn, getwt, err_fmt
from src.filesystem import get_curr_res_dpath, is_fastq
from src.fmt_readID import fmt_read_id
from src.sorter_modules.fast5 import fast5_readids, copy_read_f5_2_f5, copy_single_f5


def sort_fast5_file(f5_path, tax_annot_res_dir, sens,
        min_qual, min_qlen, logfile_path):
    """
    Function sorts FAST5 file without untwisting.

    :param f5_path: path to FAST5 file meant to be processed;
    :type f5_path: str;
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
    srt_file_dict = dict()

    trash_fpath = os.path.join(outdir_path, "qual_less_Q{}{}.fast5".format(int(min_qual),
            minlen_fmt_str))

    new_dpath = glob("{}{}*{}*".format(tax_annot_res_dir, os.sep, get_checkstr(f5_path)))[0]
    tsv_res_fpath = get_res_tsv_fpath(new_dpath)
    resfile_lines = configure_resfile_lines(tsv_res_fpath, sens)

    # File validation:
    #   RuntimeError will be raised if FAST5 file is broken.
    try:
        # File existance checking is performed while parsing CL arguments.
        # Therefore, this if-statement will trigger only if f5_path's file is not a valid HDF5 file.
        if not h5py.is_hdf5(f5_path):
            raise RuntimeError("file is not of HDF5 (i.e. not FAST5) format")
        # end if

        from_f5 = h5py.File(f5_path, 'r')

        for _ in from_f5:
            break
        # end for
    except RuntimeError as runterr:
        printl(logfile_path, err_fmt("FAST5 file is broken"))
        printl(logfile_path, "Reading the file '{}' crashed.".format(os.path.basename(fpath)))
        printl(logfile_path, "Reason: {}".format( str(runterr) ))
        printl(logfile_path, "Omitting this file...\n")
        # Return zeroes -- inc_val won't be incremented and this file will be omitted
        return (0, 0)
    # end try

    # singleFAST5 and multiFAST5 files should be processed in different ways
    # "Raw" group always in singleFAST5 root and never in multiFAST5 root
    if "Raw" in from_f5.keys():
        f5_cpy_func = copy_single_f5
    else:
        f5_cpy_func = copy_read_f5_2_f5
    # end if

    for i, read_name in enumerate(fast5_readids(from_f5)):

        try:
            hit_names, ph33_qual, q_len = resfile_lines[sys.intern(fmt_read_id(read_name))[1:]] # omit 'read_' in the beginning of FAST5 group's name
        except KeyError:
            printl(logfile_path, err_fmt("""read '{}' not found in TSV file containing taxonomic annotation.
  This TSV file: '{}'""".format(fmt_read_id(read_name), tsv_res_fpath)))
            printl(logfile_path, "Try running sorter with '-u' (--untwist-fast5') flag.\n")
            platf_depend_exit(1)
        # end try
        # If read is found in TSV file:
        q_len = SeqLength(q_len)
        if ph33_qual != '-' and ph33_qual < min_qual:
            # Get name of result FASTQ file to write this read in
            if trash_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, trash_fpath, logfile_path)
            # end if
            f5_cpy_func(from_f5, read_name, srt_file_dict[trash_fpath], logfile_path)
            seqs_fail += 1
        else:
            for hit_name in hit_names.split("&&"): # there can be multiple hits for single query sequence
                # Get name of result FASTQ file to write this read in
                sorted_file_path = os.path.join(outdir_path, "{}.fast5".format(hit_name))
                if sorted_file_path not in srt_file_dict.keys():
                    srt_file_dict = update_file_dict(srt_file_dict, sorted_file_path, logfile_path)
                # end if
                f5_cpy_func(from_f5, read_name, srt_file_dict[sorted_file_path], logfile_path)
            # end for
            seqs_pass += 1
        # end if
    # end for
    # Close all sorted files
    for file_obj in srt_file_dict.values():
        file_obj.close()
    # end for

    printl(logfile_path, "\r{} - File '{}' is sorted.".format(getwt(), os.path.basename(f5_path)))
    printn(" Working...")

    return (seqs_pass, seqs_fail)
# end def sort_fast5_file
