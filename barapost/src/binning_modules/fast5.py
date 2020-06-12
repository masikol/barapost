# -*- coding: utf-8 -*-
# This module defines functions, via which barapost-binning manipulated FAST5 data

import sys
import h5py

from src.platform import platf_depend_exit
from src.printlog import printl
from src.fmt_readID import fmt_read_id
from re import search as re_search


def fast5_readids(fast5_file):
    """
    Generator yields IDs of all reads in a FAST5 file.

    :param fast5_file: object of HDF5 FAST5 file;
    :type fast5_file: h5py.File;
    """

    if "Raw" in fast5_file.keys(): # single-FAST5 file
        yield "read_" + fmt_read_id(fast5_file.filename) # name of read is in filename
    else: # multi-FAST5 file
        for readid in fast5_file:
            if readid.startswith("read_"):
                yield readid
            # end if
        # end for
    # end if
# end def fast5_readids


def copy_read_f5_2_f5(from_f5, read_name, to_f5, logfile_path):
    """
    Function copies a read with ID 'read_name'
        from 'from_f5' multiFAST5 file to to_f5 multiFAST5 one.

    :param from_f5: FAST5 file object to copy a read from;
    :type from_f5: h5py.File;
    :param read_name: ID of a read to copy;
    :type read_name: str;
    :param to_f5: destination FAST5 file;
    :type to_f5: h5py.File;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """
    try:
        from_f5.copy(read_name, to_f5)
    except Exception as verr:
        printl(logfile_path, "\n\n ! - Error: {}".format( str(verr) ))
        printl(logfile_path, "Reason is probably the following:")
        printl(logfile_path, "  read that is copying to the result file is already in it.")
        printl(logfile_path, "ID of the read: '{}'".format(read_name))
        printl(logfile_path, "File: '{}'".format(to_f5.filename))
        platf_depend_exit(1)
    # end try
# end def copy_read_f5_2_f5


def copy_single_f5(from_f5, read_name, to_f5, logfile_path):
    """
    Function copies a read with ID 'read_name'
        from 'from_f5' singleFAST5 file to 'to_f5' multiFAST5 one.

    :param from_f5: FAST5 file object to copy a read from;
    :type from_f5: h5py.File;
    :param read_name: ID of a read to copy;
    :type read_name: str;
    :param to_f5: destination FAST5 file;
    :type to_f5: h5py.File;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """
    try:
        read_group = read_name
        to_f5.create_group(read_group) # create group in destination multi_FAST5 file

        # Copy "UniqueGlobalKey" to root of recently created group
        for ugk_subgr in from_f5["UniqueGlobalKey"]:
            from_f5.copy("UniqueGlobalKey/"+ugk_subgr, to_f5[read_group])
        # end for

        # Get data array in single-FAST5 file
        read_number_group = "Raw/Reads/"+next(iter(from_f5["Raw"]["Reads"]))
        # It's name in multi-FAST5 file
        read_number = re_search(r"(Read_[0-9]+)", read_number_group).group(1)

        # Copy group to multi-FAST5 file
        from_f5.copy(from_f5[read_number_group], to_f5[read_group])
        # Move data array to "Raw" group, as it is in multi-FAST5 files
        to_f5.move("{}/{}".format(read_group, read_number), "{}/Raw".format(read_group))

        # Copy everything else to recently created group
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

def assign_version_2(fast5_list):
    """ Assign version attribute to '2.0' -- multiFAST5"""
    for f5path in fast5_list:
        with h5py.File(f5path, 'a') as f5file:
            f5file.attrs["file_version"] = b"2.0"
        # end with
    # end for
# end def assign_version_2


def update_file_dict(srt_file_dict, new_fpath, logfile_path):
    try:
        srt_file_dict[sys.intern(new_fpath)] = h5py.File(new_fpath, 'a')
    except OSError as oserr:
        printl(logfile_path, err_fmt("error while opening one of result files"))
        printl(logfile_path, "Errorneous file: '{}'".format(new_fpath))
        printl(logfile_path,  str(oserr) )
        platf_depend_exit(1)
    # end try
    return srt_file_dict
# end def update_file_dict