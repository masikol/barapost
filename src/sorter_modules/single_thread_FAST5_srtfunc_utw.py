# -*- coding: utf-8 -*-

from src.sorter_modules.common import *

try:
    import h5py
except ImportError as imperr:
    print(err_fmt("package 'h5py' is not installed"))
    print( "Exact error description given by the interpreter: {}".format(str(imperr)) )
    print("\n  'h5py' package is necessary for FAST5 files sorting.")
    print("  Please, install it (e.g. 'pip3 install h5py').")
    print("  Tip for Linux users: you may need to install 'libhdf5-dev' with your packet manager first and then go to pip.")
    platf_depend_exit(1)
# end try


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
        printl("\n\n ! - Error: {}".format( str(verr) ))
        printl("Reason is probably the following:")
        printl("  read that is copying to the result file is already in it.")
        printl("ID of the read: '{}'".format(read_name))
        printl("File: '{}'".format(to_f5.filename))
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
        printl("\n\n ! - Error: {}".format( str(verr) ))
        printl("Reason is probably the following:")
        printl("  read that is copying to the result file is already in it.")
        printl("ID of the read: '{}'".format(read_name))
        printl("File: '{}'".format(to_f5.filename))
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

def sort_fast5_file(f5_path):
    """
    Function sorts FAST5 file with untwisting.

    :param f5_path: path to FAST5 file meant to be processed;
    :type f5_path: str;
    """

    seqs_pass = 0
    seqs_fail = 0
    srt_file_dict = dict()

    trash_fpath = os.path.join(outdir_path, "qual_less_Q{}{}.fast5".format(int(min_ph33_qual),
            minlen_fmt_str))

    from_f5 = h5py.File(f5_path, 'r') # open source FAST5
    num_reads = len(from_f5) # get number of reads in it

    # singleFAST5 and multiFAST5 files should be processed in different ways
    # "Raw" group always in singleFAST5 root and never in multiFAST5 root
    if "Raw" in from_f5.keys():
        f5_cpy_func = copy_single_f5
    else:
        f5_cpy_func = copy_read_f5_2_f5
    # end if

    try:
        readids_to_seek = list(from_f5.keys()) # list of not-sorted-yet read IDs
    except Exception as e:
        print(str(e))
        exit(0)
    # end try

    # Fill the list 'readids_to_seek'
    for read_name in fast5_readids(from_f5):
        # Get rid of "read_"
        readids_to_seek.append(sys.intern(read_name))
    # end for

    # Walk through the index
    index_f5_2_tsv = open_shelve( os.path.join(index_dirpath, index_name), 'r' )

    if not f5_path in index_f5_2_tsv.keys():
        printl(err_fmt("Source FAST5 file not found in index"))
        printl("Try to rebuild index")
        platf_depend_exit(1)
    # end if

    for tsv_path in index_f5_2_tsv[f5_path].keys():

        read_names = index_f5_2_tsv[f5_path][tsv_path]
        resfile_lines = configure_resfile_lines(tsv_path)

        for read_name in read_names:
            try:
                hit_name, ph33_qual, q_len = resfile_lines[sys.intern(fmt_read_id(read_name))]
            except KeyError:
                printl(err_fmt("missing taxonomic annotation info for read '{}'".format(fmt_read_id(read_name))))
                printl("It is stored in '{}' FAST5 file".format(f5_path))
                printl("Try to make new index file (press ENTER on corresponding prompt).")
                printl("""Or, if does not work for you, make sure that taxonomic annotation info
for this read is present in one of TSV files generated by 'prober.py' and 'barapost.py'.""")
                index_f5_2_tsv.close()
                platf_depend_exit(1)
            else:
                q_len = SeqLength(q_len)
                if q_len < min_qlen or (ph33_qual != '-' and ph33_qual < min_ph33_qual):
                    # Get name of result FASTQ file to write this read in
                    if trash_fpath not in srt_file_dict.keys():
                        srt_file_dict = update_file_dict(srt_file_dict, trash_fpath)
                    # end if
                    f5_cpy_func(from_f5, read_name, srt_file_dict[trash_fpath])
                    seqs_fail += 1
                else:
                    # Get name of result FASTQ file to write this read in
                    sorted_file_path = os.path.join(outdir_path, "{}.fast5".format(hit_name))
                    if sorted_file_path not in srt_file_dict.keys():
                        srt_file_dict = update_file_dict(srt_file_dict, sorted_file_path)
                    # end if
                    f5_cpy_func(from_f5, read_name, srt_file_dict[sorted_file_path])
                    seqs_pass += 1
                # end if
            # end try
        # end for

    index_f5_2_tsv.close()
    # Close all sorted files
    for file_obj in srt_file_dict.values():
        file_obj.close()
    # end for
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
        printl(err_fmt("error while opening one of result files"))
        printl("Errorneous file: '{}'".format(new_fpath))
        printl( str(oserr) )
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