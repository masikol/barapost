# -*- coding: utf-8 -*-
# This module defines finctions and other objects, which interact with file system.

import os
from re import search as re_search
from gzip import open as open_as_gzip
from src.platform import platf_depend_exit
from src.printlog import printl, err_fmt, getwt

# For opening plain text and gzipped files
OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)

is_gzipped = lambda file: True if file.endswith(".gz") else False

is_fastq = lambda f: True if not re_search(r".+\.f(ast)?q(\.gz)?$", f) is None else False
is_fasta = lambda f: True if not re_search(r".+\.(m)?f(ast)?a(\.gz)?$", f) is None else False
is_fast5 = lambda f: f.endswith(".fast5")


def rename_file_verbosely(file, logfile_path):
    """
    Function verbosely renames file (as well as directory) given to it.

    :param file: path to file (directory) meant to be renamed;
    :type file: str;
    :param logfile_path: path to logfile;
    :type logfile_path: str;
    """

    if not os.path.exists(file):
        return
    # end if

    # Path to "file's" parent directory
    pardir = os.path.abspath(os.path.dirname(file))
    
    # Function can rename directories
    if os.path.isdir(file):
        is_analog = lambda f: not re_search(r"{}.*(_old_[0-9]+)?$".format(os.path.basename(file)), f) is None
        word = "directory"
        name_itself = file
        ext = ""
    else:
        is_analog = lambda f: re_search(r"(.*)\..*$", os.path.basename(file)).group(1) in f
        word = "file"
        name_itself = re_search(r"(.*)\..*$", file).group(1)
        ext = re_search(r".*(\..*)$", file).group(1)
    # end if

    # Count files in 'pardir' that have analogous names as 'file' has:
    num_analog_files = len( list(filter(is_analog, os.listdir(pardir))) )

    if re_search(r"_old_[0-9]+", file) is None:
        # Append "_old_<number>"
        new_name = name_itself + "_old_" + str(num_analog_files) + ext
    else:
        # Merely substitute new number
        new_name = file.replace(re_search(r"_old_([0-9]+)", file).group(1),
            str(num_analog_files+1))
    # end if

    try:
        printl(logfile_path, '\n' + getwt() + " - Renaming old {}:".format(word))
        printl(logfile_path, "  '{}' --> '{}'".format(file, new_name))
        os.rename(file, new_name)
    except Exception as err:
        print("\n {} '{}' cannot be renamed:".format( word, str(file)) )
        print( str(err) + '\n')
        platf_depend_exit(1)
    # end try
# end def rename_file_verbosely

def remove_tmp_files(*paths):
    """
    Function removes files passed to it.
    Actually, passed arguments are paths ('str') to files meant to be removed.
    """

    for path in paths:
        if os.path.exists(path):
            try:
                os.unlink(path)
            except OSError as oserr:
                printl(err_fmt("cannot remove file '{}'").format(path))
                printl( str(oserr) )
                platf_depend_exit(1)
            # end try
        # end if
    # end for
# end def remove_tmp_files


def create_result_directory(fq_fa_path, outdir_path):
    """
    Function creates a result directory named according
        to how source FASTQ or FASTA file is named.

    :param fq_fa_path: path to source FASTQ or FASTA file;
    :type fq_fa_path: str;
    :param outdir_path: path to directory in which result_directory will be created;
    :type outdir_path: str;

    Returns 'str' path to the recently created result directory.
    """

    # dpath means "directory path"
    new_dpath = os.path.join(outdir_path, os.path.basename(fq_fa_path)) # get rid of absolute path
    new_dpath = re_search(r"(.*)\.(m)?f(ast)?(a|q)(\.gz)?$", new_dpath).group(1) # get rid of extention
    if not os.path.exists(new_dpath):
        try:
            os.makedirs(new_dpath)
        except OSError as oserr:
            printl(err_fmt("can't create result directory"))
            printl( str(oserr) )
            platf_depend_exit(1)
        # end try
    # end if
    return new_dpath
# end def create_result_directory


def get_curr_res_dpath(fq_fa_path, tax_annot_res_dir):
    """
    Function configures and returns the path to result directory for particular FASTA or FASTQ file.
    Result directory is nested in 'tax_annot_res_dir' and is named according to name of FASTA/FASTQ file.
    E.g. if file 'some_reads.fastq' is processing, it's result directory will be named 'some_reads'.

    :param fq_fa_path: path to FASTA/FASTQ file that is procesing;
    :type fq_fa_path: str;
    :param tax_annot_res_dir: path to directory with results of classification;
    :type tax_annot_res_dir: str;

    Returns path to result directory that was recenly created of 'str' type.
    """

    # dpath means "directory path"
    new_dpath = os.path.join(tax_annot_res_dir, os.path.basename(fq_fa_path)) # get rid of absolute path
    new_dpath = re_search(r"(.*)\.(m)?f(ast)?(a|q)", new_dpath).group(1) # get rid of extention

    return new_dpath
# end def get_curr_res_dir