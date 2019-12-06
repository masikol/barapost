# -*- coding: utf-8 -*-

from src.sorter_modules.common import *  # ??

from shelve import open as open_shelve

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

index_name = "fast5_to_tsvtaxann_idx"

def map_f5reads_2_taxann(fast5_list):
    """
    Function perform mapping of all reads stored in input FAST5 files
        to existing TSV files containing taxonomic annotation info.

    It creates an index DBM file.

    Generally speaking, reads from one FAST5 are spread between several
    FASTQ (and hence, TSV-taxann) files.
    Structure of our index allows to minimize times needed to read plain
    (i.e. sequential access) TSV files.
    Well, structure of our index is following:

    <DBM file>:
    {
        <path_to_FAST5_1>: {
                            <path_to_TSV_1.1>: [<read_ID_1.1.1>, <read_ID_1.1.2>, ..., <read_ID_1.1.N>],
                            <path_to_TSV_1.2>: [<read_ID_1.2.1>, <read_ID_1.2.2>, ..., <read_ID_1.2.N>],
                            ...
                            <path_to_TSV_1.M>: [<read_ID_1.M.1>, <read_ID_1.M.2>, ..., <read_ID_1.M.N>]
                         },
        <path_to_FAST5_2>: {
                            <path_to_TSV_1.1>: [<read_ID_1.1.1>, <read_ID_1.1.2>, ..., <read_ID_1.1.N>],
                            <path_to_TSV_1.2>: [<read_ID_1.2.1>, <read_ID_1.2.2>, ..., <read_ID_1.2.N>],
                            ...
                            <path_to_TSV_2.M>: [<read_ID_2.M.1>, <read_ID_2.M.2>, ..., <read_ID_2.M.N>]
                         },
        ...
        <path_to_FAST5_K>: {
                            <path_to_TSV_K.1>: [<read_ID_K.1.1>, <read_ID_K.1.2>, ..., <read_ID_K.1.N>],
                            <path_to_TSV_K.2>: [<read_ID_K.2.1>, <read_ID_K.2.2>, ..., <read_ID_K.2.N>],
                            ...
                            <path_to_TSV_K.M>: [<read_ID_K.M.1>, <read_ID_K.M.2>, ..., <read_ID_K.M.N>]
                         },
    }
    """

    # Get all directories nested in 'tax_annot_res_dir'
    taxann_dir_lst = list(filter(lambda f: True if os.path.isdir(f) else False,
        glob( os.path.join(tax_annot_res_dir, "*") )))

    # Exclude "local_database" and "fast5_to_fastq_idx" from this list
    for dir_to_exclude in (index_name, "local_database"):
        ldb_dir_path = os.path.join(tax_annot_res_dir, dir_to_exclude)
        if ldb_dir_path in taxann_dir_lst:
            taxann_dir_lst.remove(ldb_dir_path)
        # end if
    # end for

    # Get path to TSV files containing taxonomy annotation info
    tsv_taxann_lst = list()
    for taxann_dir in taxann_dir_lst:
        putative_tsvs = glob("{}{}*.tsv".format(taxann_dir, os.sep))
        if len(putative_tsvs) == 1:
            tsv_taxann_lst.append(putative_tsvs[0])
        elif len(putative_tsvs) == 0:
            printl("""Warning!  There is no taxonomic annotation info in the following directory:
'{}'
Omitting this directory.\n""".format(taxann_dir))
        else:
            printl("""Error!  Multiple TSV files in the following directory:
'{}'
Please, remove extra files and leave only one, which contains actual taxononic annotation info.""".format(taxann_dir))
            platf_depend_exit(1)
        # end if
    # end for
    del taxann_dir_lst

    printl("{} - Untwisting started.".format(getwt()))

    # Open index files overwriting existing data ('n' parameter)
    index_f5_2_tsv = open_shelve( os.path.join(index_dirpath, index_name), 'n' )

    global inc_val
    inc_val = 0
    get_inc_val = lambda: inc_val # merely return this value (1 thread)

    # Launch printer
    printer = Thread(target=status_printer, args=(get_inc_val, stop)) # create thread
    stop = Event()
    stop.set() # raise the flag
    printer.start() # start waiting

    # Iterate over FAST5 files
    for j, f5_path in enumerate(fast5_list):

        f5_file = h5py.File(f5_path, 'r')# open FAST5 file
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
                    if fmt_read_id(readid) in readids_in_tsv:
                        # If not first -- write data to dict (and to index later)
                        try:
                            idx_dict[tsv_taxann_fpath].append(readid) # append to existing list
                        except KeyError:
                            idx_dict[tsv_taxann_fpath] = [readid] # create a new list
                        finally:
                            readids_to_seek.remove(readid)
                            inc_val += 1
                        # end try
                    # end if
                # end for
            # end with
            if len(readids_to_seek) == 0:
                break
            # end if
        # end for

        # If after all TSV is checked but nothing have changed -- we miss taxonomic annotation
        #     for some reads! And we will write their IDs to 'missing_reads_lst.txt' file.
        if len(readids_to_seek) == len_before:
            printl(err_fmt("reads from FAST5 file not found"))
            printl("FAST5 file: '{}'".format(f5_path))
            printl("Some reads reads have not undergone taxonomic annotation.")
            missing_log = "missing_reads_lst.txt"
            printl("List of missing reads are in following file:\n  '{}'\n".format(missing_log))
            with open(missing_log, 'w') as missing_logfile:
                missing_logfile.write("Missing reads from file '{}':\n\n".format(f5_path))
                for readid in readids_to_seek:
                    missing_logfile.write(fmt_read_id(readid) + '\n')
                # end for
            index_f5_2_tsv.close()
            try:
                for path in glob( os.path.join(index_dirpath, '*') ):
                    os.unlink(path)
                # end for
                os.rmdir(index_dirpath)
            except OSError as oserr:
                printl("error while removing index directory: {}".format(oserr))
            finally:
                platf_depend_exit(3)
            # end try
        # end if

        # Update index
        index_f5_2_tsv[f5_path] = idx_dict
    # end for

    # Stop printer
    stop = True # lower the flag
    printer.join()
    printl()

    index_f5_2_tsv.close()
    printl("{} - Untwisting is completed.".format(getwt()))
    printl('-'*20+'\n')
# end def map_f5reads_2_taxann