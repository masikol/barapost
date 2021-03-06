# -*- coding: utf-8 -*-
# This module defines functions, via which barapost-binning manipulated FAST5 data

import re
import sys
import h5py

from src.platform import platf_depend_exit
from src.printlog import printlog_error, printlog_error_time
from src.fmt_read_id import fmt_read_id


def fast5_readids(fast5_file):
    # Generator yields IDs of all reads in a FAST5 file.
    # :param fast5_file: object of HDF5 FAST5 file;
    # :type fast5_file: h5py.File;

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


def copy_read_f5_2_f5(from_f5, read_name, to_f5):
    # Function copies a read with ID 'read_name'
    #     from 'from_f5' multiFAST5 file to to_f5 multiFAST5 one.
    #
    # :param from_f5: FAST5 file object to copy a read from;
    # :type from_f5: h5py.File;
    # :param read_name: ID of a read to copy;
    # :type read_name: str;
    # :param to_f5: destination FAST5 file;
    # :type to_f5: h5py.File;

    if not to_f5 is None: # handle no_trash
        try:
            from_f5.copy(read_name, to_f5)
        except ValueError as err:
            printlog_error_time("Error: `{}`".format( str(err) ))
            printlog_error("Reason is probably the following:")
            printlog_error("  read that is copying to the result file is already in this file.")
            printlog_error("ID of the read: `{}`".format(read_name))
            printlog_error("File: `{}`".format(to_f5.filename))
            return
        # end try
    # end if
# end def copy_read_f5_2_f5


def copy_single_f5(from_f5, read_name, to_f5):
    # Function copies a read with ID 'read_name'
    #     from 'from_f5' singleFAST5 file to 'to_f5' multiFAST5 one.
    #
    # :param from_f5: FAST5 file object to copy a read from;
    # :type from_f5: h5py.File;
    # :param read_name: ID of a read to copy;
    # :type read_name: str;
    # :param to_f5: destination FAST5 file;
    # :type to_f5: h5py.File;

    # Handle no_trash
    if to_f5 is None:
        return
    # end if

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
        read_number = re.search(r"(Read_[0-9]+)", read_number_group).group(1)

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
    except ValueError as err:
        printlog_error_time("Error: `{}`".format( str(err) ))
        printlog_error("Reason is probably the following:")
        printlog_error("  read that is copying to the result file is already in this file.")
        printlog_error("ID of the read: `{}`".format(read_name))
        printlog_error("File: `{}`".format(to_f5.filename))
        return
    # end try
# end def copy_single_f5


def assign_version_2(fast5_list):
    # Assign version attribute to '2.0' -- multiFAST5
    for f5path in fast5_list:
        with h5py.File(f5path, 'a') as f5file:
            f5file.attrs["file_version"] = b"2.0"
        # end with
    # end for
# end def assign_version_2


def update_file_dict(srt_file_dict, new_fpath):
    try:
        if not new_fpath is None:
            srt_file_dict[sys.intern(new_fpath)] = h5py.File(new_fpath, 'a')
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
