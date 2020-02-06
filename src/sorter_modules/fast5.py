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


from src.printlog import printl
from re import search as re_search


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


def copy_single_f5(from_f5, read_name, to_f5, logfile_path):
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