# -*- coding: utf-8 -*-
# Module defines function that creates index mapping reads in FAST5 files to their
#   taxonomic annotation in single thread.

import os
import h5py

from glob import glob
from shelve import open as open_shelve

from src.fmt_read_id import fmt_read_id
from src.platform import platf_depend_exit
from src.binning_modules.fast5 import fast5_readids
from src.printlog import printlog_error, printlog_error_time

index_name = "fast5_to_tsvtaxann_idx"


#   Structure of index:
# Generally speaking, reads from one FAST5 are spread between several
# FASTQ (and hence, TSV-taxann) files.
# Structure of our index allows to minimize times needed to read plain
# (i.e. sequential access) TSV files.
# Well, structure of our index is following:

# <DBM file>:
# {
#     <path_to_FAST5_1>: {
#        <path_to_TSV_1.1>: [<read_ID_1.1.1>, <read_ID_1.1.2>, ..., <read_ID_1.1.N>],
#        <path_to_TSV_1.2>: [<read_ID_1.2.1>, <read_ID_1.2.2>, ..., <read_ID_1.2.N>],
#        ...
#        <path_to_TSV_1.M>: [<read_ID_1.M.1>, <read_ID_1.M.2>, ..., <read_ID_1.M.N>]
#        },
#     <path_to_FAST5_2>: {
#        <path_to_TSV_1.1>: [<read_ID_1.1.1>, <read_ID_1.1.2>, ..., <read_ID_1.1.N>],
#        <path_to_TSV_1.2>: [<read_ID_1.2.1>, <read_ID_1.2.2>, ..., <read_ID_1.2.N>],
#        ...
#        <path_to_TSV_2.M>: [<read_ID_2.M.1>, <read_ID_2.M.2>, ..., <read_ID_2.M.N>]
#        },
#     ...
#     <path_to_FAST5_K>: {
#        <path_to_TSV_K.1>: [<read_ID_K.1.1>, <read_ID_K.1.2>, ..., <read_ID_K.1.N>],
#        <path_to_TSV_K.2>: [<read_ID_K.2.1>, <read_ID_K.2.2>, ..., <read_ID_K.2.N>],
#        ...
#        <path_to_TSV_K.M>: [<read_ID_K.M.1>, <read_ID_K.M.2>, ..., <read_ID_K.M.N>]
#        }
# }


def map_f5reads_2_taxann(f5_path, tsv_taxann_lst, tax_annot_res_dir):
    # Function perform mapping of all reads stored in input FAST5 files
    #     to existing TSV files containing taxonomic annotation info.
    #
    # It creates an DBM index file.
    #
    # :param f5_path: path to current FAST5 file;
    # :type f5_path: str;
    # :param tsv_taxann_lst: list of path to TSV files that contain taxonomic annotation;
    # :type tsv_taxann_lst: list<str>;
    # :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    # :type tax_annot_res_dir: str;

    index_dirpath = os.path.join(tax_annot_res_dir, index_name) # name of directory that will contain indicies

    # File validation:
    #   RuntimeError will be raised if FAST5 file is broken.
    try:
        # File existance checking is performed while parsing CL arguments.
        # Therefore, this if-statement will trigger only if f5_path's file is not a valid HDF5 file.
        if not h5py.is_hdf5(f5_path):
            raise RuntimeError("file is not of HDF5 (i.e. not FAST5) format")
        # end if

        f5_file = h5py.File(f5_path, 'r')

        for _ in f5_file:
            break
        # end for
    except RuntimeError as runterr:
        printlog_error_time("FAST5 file is broken")
        printlog_error("Reading the file `{}` crashed.".format(os.path.basename(f5_path)))
        printlog_error("Reason: {}".format( str(runterr) ))
        printlog_error("Omitting this file...")
        print()
        return
    # end try

    readids_to_seek = list(fast5_readids(f5_file))
    idx_dict = dict() # dictionary for index

    # This saving is needed to compare with 'len(readids_to_seek)'
    #    after all TSV will be looked through in order to
    #    determine if some reads miss taxonomic annotation.
    len_before = len(readids_to_seek)

    # Iterate over TSV-taaxnn file
    for tsv_taxann_fpath in tsv_taxann_lst:

        with open(tsv_taxann_fpath, 'r') as taxann_file:

            # Get all read IDs in current TSV
            readids_in_tsv = list( map(lambda l: l.split('\t')[0], taxann_file.readlines()) )

            # Iterate over all other reads in current FAST5
            #    ('reversed' is necessary because we remove items from list in this loop)
            for readid in reversed(readids_to_seek):
                fmt_id = fmt_read_id(readid)[1:]
                if fmt_id in readids_in_tsv:
                    # If not first -- write data to dict (and to index later)
                    try:
                        idx_dict[tsv_taxann_fpath].append("read_"+fmt_id) # append to existing list
                    except KeyError:
                        idx_dict[tsv_taxann_fpath] = ["read_" + fmt_id] # create a new list
                    finally:
                        readids_to_seek.remove(readid)
                    # end try
                # end if
            # end for
        # end with
        if len(readids_to_seek) == 0:
            break
        # end if
    # end for

    # Save info about reads, for which classification if not found
    #   in any of classification files
    if len(readids_to_seek) != 0:
        not_fount_key = 'CLASSIF_NOT_FOUND'
        idx_dict[not_fount_key] = list()
        for readid in readids_to_seek:
            fmt_id = fmt_read_id(readid)[1:]
            idx_dict[not_fount_key].append("read_" + fmt_id)
        # end for
    # end if

    # If after all TSV is checked but nothing have changed -- we miss taxonomic annotation
    #     for some reads! And we will write their IDs to 'missing_reads_lst.txt' file.
    if len(readids_to_seek) == len_before:
        printlog_error_time("reads from FAST5 file not found")
        printlog_error("FAST5 file: `{}`".format(f5_path))
        printlog_error("Some reads have not undergone taxonomic annotation.")
        missing_log = "missing_reads_lst.txt"
        printlog_error("List of missing reads are in following file:")
        printlog_error("{}".format(missing_log))
        with open(missing_log, 'w') as missing_logfile:
            missing_logfile.write("Missing reads from file '{}':\n\n".format(f5_path))
            for readid in readids_to_seek:
                missing_logfile.write(fmt_read_id(readid) + '\n')
            # end for
    # end if

    try:
        # Open index files appending to existing data ('c' parameter)
        with open_shelve( os.path.join(index_dirpath, index_name), 'c' ) as index_f5_2_tsv:
            # Update index
            index_f5_2_tsv[f5_path] = idx_dict
        # end with
    except OSError as oserr:
        printlog_error_time("Error: cannot create index file `{}`"\
            .format(os.path.join(index_dirpath, index_name)))
        printlog_error( str(oserr) )
        platf_depend_exit(1)
    # end try
# end def map_f5reads_2_taxann
