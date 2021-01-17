# -*- coding: utf-8 -*-
# Module defines functions necessary for binning FAST5 files "directly": without index.

import h5py

import os
import sys
from glob import glob
import logging

from src.binning_modules.binning_spec import get_checkstr, get_res_tsv_fpath, configure_resfile_lines
from src.binning_modules.fast5 import update_file_dict
from src.binning_modules.fast5 import fast5_readids, copy_read_f5_2_f5, copy_single_f5

from src.binning_modules.filters import get_QL_filter, get_QL_trash_fpath
from src.binning_modules.filters import get_align_filter, get_align_trash_fpath

from src.platform import platf_depend_exit
from src.printlog import printn, printlog_info, printlog_error, printlog_error_time, printlog_info_time
from src.fmt_readID import fmt_read_id


def bin_fast5_file(f5_path, tax_annot_res_dir, sens, min_qual, min_qlen,
    min_pident, min_coverage, no_trash):
    # Function bins FAST5 file without untwisting.
    #
    # :param f5_path: path to FAST5 file meant to be processed;
    # :type f5_path: str;
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

    srt_file_dict = dict()

    new_dpath = glob("{}{}*{}*".format(tax_annot_res_dir, os.sep, get_checkstr(f5_path)))[0]
    tsv_res_fpath = get_res_tsv_fpath(new_dpath)
    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy", "taxonomy.tsv")
    resfile_lines = configure_resfile_lines(tsv_res_fpath, sens, taxonomy_path)

    # Make filter for quality and length
    QL_filter = get_QL_filter(f5_path, min_qual, min_qlen)
    # Configure path to trash file
    if not no_trash:
        QL_trash_fpath = get_QL_trash_fpath(f5_path, outdir_path, min_qual, min_qlen,)
    else:
        QL_trash_fpath = None
    # end if

    # Make filter for identity and coverage
    align_filter = get_align_filter(min_pident, min_coverage)
    # Configure path to this trash file
    if not no_trash:
        align_trash_fpath = get_align_trash_fpath(f5_path, outdir_path, min_pident, min_coverage)
    else:
        align_trash_fpath = None
    # end if

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
        printlog_error_time("FAST5 file is broken")
        printlog_error("Reading the file `{}` crashed.".format(os.path.basename(f5_path)))
        printlog_error("Reason: {}".format( str(runterr) ))
        printlog_error("Omitting this file...")
        print()
        # Return zeroes -- inc_val won't be incremented and this file will be omitted
        return (0, 0, 0)
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
            hit_names, *vals_to_filter = resfile_lines[sys.intern(fmt_read_id(read_name))[1:]] # omit 'read_' in the beginning of FAST5 group's name
        except KeyError:
            printlog_error_time("Error: read `{}` not found in TSV file containing taxonomic annotation.")
            printlog_error("This TSV file: `{}`".format(fmt_read_id(read_name), tsv_res_fpath))
            printlog_error("Try running barapost-binning with `-u` (`--untwist-fast5`) flag.\n")
            platf_depend_exit(1)
        # end try
        # If read is found in TSV file:
        if not QL_filter(vals_to_filter):
            QL_seqs_fail += 1
            # Get name of result FASTQ file to write this read in
            if QL_trash_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, QL_trash_fpath)
            # end if
            f5_cpy_func(from_f5, read_name, srt_file_dict[QL_trash_fpath])
        elif not align_filter(vals_to_filter):
            align_seqs_fail += 1
            # Get name of result FASTQ file to write this read in
            if QL_trash_fpath not in srt_file_dict.keys():
                srt_file_dict = update_file_dict(srt_file_dict, align_trash_fpath)
            # end if
            f5_cpy_func(from_f5, read_name, srt_file_dict[align_trash_fpath])
        else:
            for hit_name in hit_names.split("&&"): # there can be multiple hits for single query sequence
                # Get name of result FASTQ file to write this read in
                binned_file_path = os.path.join(outdir_path, "{}.fast5".format(hit_name))
                if binned_file_path not in srt_file_dict.keys():
                    srt_file_dict = update_file_dict(srt_file_dict, binned_file_path)
                # end if
                f5_cpy_func(from_f5, read_name, srt_file_dict[binned_file_path])
            # end for
            seqs_pass += 1
        # end if
    # end for

    from_f5.close()

    # Close all binned files
    for file_obj in filter(lambda x: not x is None, srt_file_dict.values()):
        file_obj.close()
    # end for


    sys.stdout.write('\r')
    printlog_info_time("File `{}` is binned.".format(os.path.basename(f5_path)))
    printn(" Working...")

    return (seqs_pass, QL_seqs_fail, align_seqs_fail)
# end def bin_fast5_file
