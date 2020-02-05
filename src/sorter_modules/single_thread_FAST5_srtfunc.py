# -*- coding: utf-8 -*-

from src.platform import platf_depend_exit

try:
    import h5py
except ImportError as imperr:
    print("Package 'h5py' is not installed")
    print( "Exact error description given by the interpreter: {}".format(str(imperr)) )
    print("\n  'h5py' package is necessary for FAST5 files sorting.")
    print("  Please, install it (e.g. 'pip3 install h5py').")
    print("  Tip for Linux users: you may need to install 'libhdf5-dev' with your packet manager first and then go to pip.")
    platf_depend_exit(1)
# end try

import sys
from glob import glob

from src.sorter_modules.sorter_spec import *

from src.printlog import printl, printn, getwt, err_fmt
from src.filesystem import get_curr_res_dpath, is_fastq
from src.fmt_readID import fmt_read_id


def copy_read_f5_2_f5(from_f5, read_name, to_f5):
    """
    Function copies a read with ID 'read_name'
        from 'from_f5' multiFAST5 file to to_f5 multiFAST5 one.

    :param from_f5: FAST5 file object to copy a read from;
    :type from_f5: h5py.File;
    :param read_name: ID of a read to copy;
    :type read_name: str;
    :param to_f5: destination FAST5 file;
    :type to_f5: h5py.File;
    """
    try:
        from_f5.copy(read_name, to_f5)

    except ValueError as verr:
        printl(logfile_path, "\n\n ! - Error: {}".format( str(verr) ))
        printl(logfile_path, "Reason is probably the following:")
        printl(logfile_path, "  read that is copying to the result file is already in it.")
        printl(logfile_path, "ID of the read: '{}'".format(read_name))
        printl(logfile_path, "File: '{}'".format(to_f5.filename))
        platf_depend_exit(1)
    # end try
# end def copy_read_f5_2_f5


def copy_single_f5(from_f5, read_name, to_f5):
    """
    Function copies a read with ID 'read_name'
        from 'from_f5' singleFAST5 file to to_f5 multiFAST5 one.

    :param from_f5: FAST5 file object to copy a read from;
    :type from_f5: h5py.File;
    :param read_name: ID of a read to copy;
    :type read_name: str;
    :param to_f5: destination FAST5 file;
    :type to_f5: h5py.File;
    """
    try:
        read_group = read_name
        to_f5.create_group(read_group)

        for ugk_subgr in from_f5["UniqueGlobalKey"]:
                from_f5.copy("UniqueGlobalKey/"+ugk_subgr, to_f5[read_group])
        # end for

        read_number_group = "Raw/Reads/"+next(iter(from_f5["Raw"]["Reads"]))
        read_number = re_search(r"(Read_[0-9]+)", read_number_group).group(1)
        from_f5.copy(from_f5[read_number_group], to_f5[read_group])
        to_f5.move("{}/{}".format(read_group, read_number), "{}/Raw".format(read_group))

        for group in from_f5:
            if group != "Raw" and group != "UniqueGlobalKey":
                from_f5.copy(group, to_f5["/{}".format(read_group)])
            # end if
        # end for
    except ValueError as verr:
        printl(logfile_path, "\n\n ! - Error: {}".format( str(verr) ))
        printl(logfile_path, "Reason is probably the following:")
        printl(logfile_path, "  read that is copying to the result file is already in it.")
        printl(logfile_path, "ID of the read: '{}'".format(read_name))
        printl(logfile_path, "File: '{}'".format(to_f5.filename))
        platf_depend_exit(1)
    # end try
# end def copy_single_f5


def fast5_readids(fast5_file):

    if "Raw" in fast5_file.keys():
        yield "read_" + fmt_read_id(fast5_file.filename)
        return
    else:
        for readid in fast5_file:
            if readid.startswith("read_"):
                yield readid
            # end if
        # end for
    # end if
    return
# end def fast5_readids


def sort_fast5_file(f5_path, tax_annot_res_dir, sens,
        min_qual, min_qlen, logfile_path):
    """
    Function sorts FAST5 file without untwisting.

    :param f5_path: path to FAST5 file meant to be processed;
    :type f5_path: str;
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
                srt_file_dict = update_file_dict(srt_file_dict, trash_fpath)
            # end if
            f5_cpy_func(from_f5, read_name, srt_file_dict[trash_fpath])
            seqs_fail += 1
        else:
            for hit_name in hit_names.split("&&"):
                # Get name of result FASTQ file to write this read in
                sorted_file_path = os.path.join(outdir_path, "{}.fast5".format(hit_name))
                if sorted_file_path not in srt_file_dict.keys():
                    srt_file_dict = update_file_dict(srt_file_dict, sorted_file_path)
                # end if
                f5_cpy_func(from_f5, read_name, srt_file_dict[sorted_file_path])
            # end for
            seqs_pass += 1
        # end if
    # end for
    # Close all sorted files
    for file_obj in srt_file_dict.values():
        file_obj.close()
    # end for

    printl(logfile_path, "\rFile '{}' is sorted.".format(os.path.basename(f5_path)))
    printn(" Working...")

    return (seqs_pass, seqs_fail)
# end def sort_fast5_file


def update_file_dict(srt_file_dict, new_fpath):
    try:
        if new_fpath.endswith(".fast5"):
            srt_file_dict[sys.intern(new_fpath)] = h5py.File(new_fpath, 'a')
        else:
            srt_file_dict[sys.intern(new_fpath)] = open(new_fpath, 'a')
        # end if
    except OSError as oserr:
        printl(logfile_path, err_fmt("error while opening one of result files"))
        printl(logfile_path, "Errorneous file: '{}'".format(new_fpath))
        printl(logfile_path,  str(oserr) )
        platf_depend_exit(1)
    # end try
    return srt_file_dict
# end def update_file_dict


def assign_version_2(fast5_list):
    # Assign version attribute to '2.0' -- multiFAST5
    for f5path in fast5_list:
        with h5py.File(f5path, 'a') as f5file:
            f5file.attrs["file_version"] = b"2.0"
        # end with
    # end for
# end def assign_version_2