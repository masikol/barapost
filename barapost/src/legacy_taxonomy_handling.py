# -*- coding: utf-8 -*-
# This module defines functions, which check if legacy taxonomy file is used
#   and reformats it to the new format.

import os

import src.taxonomy as taxonomy
from src.printlog import printlog_info, printlog_error, printn, log_info
from src.platform import platf_depend_exit

def check_deprecated_taxonomy(classif_dir):

    legacy_mode = False
    legacy_tax_path = os.path.join(classif_dir, "taxonomy", "taxonomy")

    if not os.path.exists(legacy_tax_path):
        pass # continue silently
    else:
        print()
        printlog_info("Legacy taxonomy file detected: `{}`.".format(legacy_tax_path))
        printlog_info("It will be reformatted to new format -- to plain TSV.")

        _reformat_legacy_file(legacy_tax_path)
    # end if
# end def check_deprecated_taxonomy


def _reformat_legacy_file(legacy_tax_path):

    import shelve

    # Check if this file is corrupted
    try:
        with shelve.open(legacy_tax_path, 'r') as tax_file:
            pass
        # end with
    except Exception as err:
        printlog_error("Legacy taxonomy file appears to be corrupted.")
        printlog_error("This error might be fatal.")
        str_err = str(err)
        if "dbm.gnu" in str_err and "module is not" in str_err:
            printlog_error("Installing `python3-gdbm` might solve this problem.")
        else:
            printlog_error("The program can't recover taxonomy from the broken file.")
            printlog_error("Seems, you have to annotate your sequences again.")
            printlog_error("Sorry for that :(")
        # end if
        platf_depend_exit(1)
    # end try

    new_tax_path = "{}.tsv".format(legacy_tax_path)

    taxonomy.init_tax_file(new_tax_path)

    printn("Reformatting: `{}` ->".format(legacy_tax_path))
    log_info("Reformatting: `{}` ->".format(legacy_tax_path))

    with shelve.open(legacy_tax_path, 'r') as old_tax_file, open(new_tax_path, 'w') as new_tax_file:
        for acc, taxonomy in old_tax_file.items():
            if isinstance(taxonomy, tuple):
                tax_str = taxonomy.config_taxonomy_str(taxonomy)
                new_tax_file.write("{}\n".format('\t'.join( (acc, tax_str) )))
            elif isinstance(taxonomy, str):
                new_tax_file.write("{}\n".format('\t'.join( (acc, taxonomy) )))
            else:
                # Execution must not reach here
                printlog_error("\nFatal error 8755.")
                printlog_error("Please, contact the developer -- it is his fault.")
                platf_depend_exit(8755)
            # end if
        # end for
    # end with

    printlog_info(" `<same_dir>/{}`".format(os.path.basename(new_tax_path)))

    try:
        renamed_legacy_file = "{}_deprecated".format(legacy_tax_path)
        os.rename(legacy_tax_path, renamed_legacy_file)
    except OSError as err:
        printlog_error("Cannot rename legacy taxonomy file `{}`:".format(legacy_tax_path))
        printlog_error(str(err))
        printlog_error("But it's not a problem -- we will proceed with our work.")
    else:
        printlog_info("Renamed: `{}` -> `<same_dir>/{}`".format(legacy_tax_path,
            os.path.basename(renamed_legacy_file)))
    # end try

    printlog_info("Legacy taxonomy file is reformatted to TSV format.")
# end def _reformat_legacy_file
