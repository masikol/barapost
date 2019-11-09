#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2019.11.10"

# |===== Check python interpreter version =====|

from sys import version_info as verinf

if verinf.major < 3:
    print( "\nYour python interpreter version is " + "%d.%d" % (verinf.major, verinf.minor) )
    print("   Please, use Python 3!\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    raw_input("Press ENTER to exit:")
    exit(1)
# end if

from sys import platform

def platf_depend_exit(exit_code):
    """
    A function that asks to press ENTER on Windows
        and exits.

    :type exit_code: int;
    """
    if platform.startswith("win"):
        input("Press ENTER to exit:")
    # end if
    exit(exit_code)
# end def platf_depend_exit

help_msg = """\nBCsummarizer.py
  Version {}; {} edition;\n
DESCRIPTION:\n
BCsummarizer.py -- this script is designed for generating a brief summary of basecalling.
  It determines, in which FASTQ files are reads from FAST5 files placed.\n
This script can be useful, because basecallers (popular Guppy, in particular)
  often missasign names of input FAST5 and output FASTQ files. In result, source FAST5
  and basecalled FASTQ files contain different reads although their names match one another.\n
Therefore, straitforward sorting of FAST5 files, that relies on names of "corresponding"
  FASTQ files (that have ondergone taxonomic annotation) is, in general, impossible.
----------------------------------------------------------\n
OPTIONS:\n
    -h (--help) --- show help message;\n
    -v (--version) --- show version;\n
    -5 (--fast5-dir) --- directory that contains FAST5 files
        meant to be processed. It may contain not only FAST5 files;\n
    -q (--fastq-dir) --- directory that contains FASTQ files
        meant to be processed. It may contain not only FASTQ files.
        FASTQ files can be gzipped;\n
    -o (--outfile) --- output summary file;
----------------------------------------------------------\n
EXAMPLES:\n
1. FAST5 files are in directory 'F5_dir'. Basecalled FASTQ files
  are in directory 'FQ_dir':\n
  ./BCsummarizer.py -5 F5_dir -q FQ_dir\n
2. FAST5 and basecalled FASTQ files are in the working directory.
  Write results in the file '/tmp/seq_summ.txt':\n
  ./BCsummarizer -5 ./ -q ./ -o /tmp/seq_summ.txt
""".format(__version__, __last_update_date__)

from sys import argv

# First search for information-providing options:

if "-h" in argv[1:] or "--help" in argv[1:]:
    print(help_msg)
    platf_depend_exit(0)
# end if

if "-v" in argv[1:] or "--version" in argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

def err_fmt(text):
    """Function for configuring error messages"""
    return "\n   \a!! - ERROR: " + text + '\n'
# end def print_error

try:
    import h5py
except ImportError as imperr:
    print(err_fmt("package 'h5py' is not installed"))
    print( "Exact error description given by the interpreter: {}".format(str(imperr)) )
    print("\n  'h5py' package is necessary for working with FAST5 files.")
    print("  Please, install it (e.g. 'pip3 install h5py').")
    print("  Tip for Linux users: you may need to install 'libhdf5-dev' with your packet manager first and then go to pip.")
    platf_depend_exit(1)
# end try


# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, sleep, gmtime
start_time = time()


def get_work_time():
    return strftime("%H:%M:%S", gmtime( time() - start_time ))
# end def get_work_time


# Get start time
from datetime import datetime
now = datetime.now().strftime("%Y-%m-%d %H.%M.%S")
# -------------------

from re import search as re_search
from glob import glob
import os
from sys import intern
from getopt import gnu_getopt, GetoptError


try:
    opts, args = gnu_getopt(argv[1:], "hv5:q:o:", ["help", "version", "fast5-dir", "fastq-dir", "outfile"])
except GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

fast5_dir = None
fastq_dir = None
outfpath = "basecall_summary.txt"

is_fqgz = lambda f: False if re_search(r".*\.f(ast)?q(\.gz)$", f) is None else True

for opt, arg in opts:

    if opt in ("-5", "--fast5-dir"):
        if not os.path.isdir(arg):
            print(err_fmt("'{}' is not a directory".format(arg)))
            platf_depend_exit(1)
        # end if
        fast5_dir = arg

        fast5_lst = glob(os.path.join(fast5_dir, "*.fast5"))
        if len(fast5_lst) == 0:
            print(err_fmt("there are no FAST5 files in '{}' directory".format(fast5_dir)))
        # end if

        fast5_lst.sort()
        fast5_lst = list(map(os.path.abspath, fast5_lst))
    # end if

    if opt in ("-q", "--fastq-dir"):
        if not os.path.isdir(arg):
            print(err_fmt("'{}' is not a directory".format(arg)))
            platf_depend_exit(1)
        # end if
        fastq_dir = arg

        fastq_lst = list(filter(is_fqgz, glob(os.path.join(fastq_dir, "*"))))
        if len(fastq_lst) == 0:
            print(err_fmt("there are no FASTQ files in '{}' directory".format(fastq_dir)))
        # end if

        fastq_lst.sort()
        fastq_lst = list(map(os.path.abspath, fastq_lst))

        # These paths (to FASTQ files) will be used as dict keys, therefore intern them
        fastq_lst = list(map(intern, fastq_lst))
    # end if

    if opt in ("-o", "--outfile"):
        outfpath = arg
        try:
            with open(outfile, 'w') as outfile:
                pass
            # end with
        except OSError as oserrr:
            print(err_fmt("cannot create outfile"))
            print( str(oserr) )
        # end try
    # end if
# end for

if fast5_dir is None:
    print(err_fmt("directory that contains FAST5 files ('-5' option) is not specified"))
    platf_depend_exit(1)
# end if

if fastq_dir is None:
    print(err_fmt("directory that contains FASTQ files ('-q' option) is not specified"))
    platf_depend_exit(1)
# end if

from sys import stdout as sys_stdout
def printn(text):
    """
    Function prints text to the console without adding '\\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
    """
    sys_stdout.write(text)
    sys_stdout.flush()
# end def printn

print(("\n |=== BCsummarizer.py (version {}) ===|\n".format(__version__)))

print( get_work_time() + " ({}) ".format(strftime("%Y.%m.%d %H:%M:%S", localtime(start_time))) + "- Start working\n")

printn("{} - Primary validation...".format(get_work_time()))
for fpath in fast5_lst:
    try:
        # File existance checking is performed while parsing CL arguments.
        # Therefore, this if-statement will trigger only if f5_path's file is not a valid HDF5 file.
        if not h5py.is_hdf5(fpath):
            raise RuntimeError("file is not of HDF5 (i.e. not FAST5) format")
        # end if
        # Open the file

        f5_file = h5py.File(fpath, 'r')
        # Validation of the file:
        #   RuntimeError will be raised if FAST5 file is broken.
        for _ in f5_file:
            break
        # end for
    except OSError as oserr:
        print(err_fmt("cannot open FAST5 file"))
        print("Opening file '{}' crashed:".format(os.path.basename(fpath)))
        print("Reason: {}".format( str(oserr) ))
        platf_depend_exit(1)
    except RuntimeError as runterr:
        print(err_fmt("FAST5 file is broken"))
        print("Reading of the file '{}' crashed.".format(os.path.basename(fpath)))
        print("Reason: {}".format( str(runterr) ))
        platf_depend_exit(5)
    finally:
        f5_file.close()
    # end try
# end for
print("\r{} - Primary validation... {}\n".format(get_work_time(), "Ok"))

print("""Results will be written to the following file:
  '{}'\n""".format(outfpath))

from gzip import open as open_as_gzip
is_gzipped = lambda f: True if f.endswith(".gz") else False

OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line,   # format text line
    lambda line: line.decode("utf-8")  # format gzipped line
)


def get_fastq_readids(fastq_fpath):

    how_to_open = OPEN_FUNCS[ is_gzipped(fastq_fpath) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq_fpath) ]

    with how_to_open(fastq_fpath) as fastq_file:
        line = fmt_func(fastq_file.readline())
        yield line.partition(' ')[0].replace('@', '')
        i = 1
        while line != "":
            line = fmt_func(fastq_file.readline())
            i += 1
            if i == 4:
                i = 0
            elif i == 1:
                yield line.partition(' ')[0].replace('@', '')
            # end if
        # end while
    # end with
# end def get_fastq_readids


final_verdict = "\nAll reads from all FAST5 files are found"
total_nlost_reads = 0

outfile = open(outfpath, 'w')
for fast5_fpath in fast5_lst:
    print("{} - '{}': start processing".format(get_work_time(), os.path.basename(fast5_fpath)))

    f5_file = h5py.File(fast5_fpath, 'r')
    num_reads = len(f5_file)

    readids_to_seek = list() # here all read IDs will be stored

    # Store add read IDs in a list
    for read_name in f5_file:
        if not read_name.startswith("read_"):
            printl("\nName of read '{}' from FAST5 file has format that is unforseen by the developer.".format(read_name))
            printl("Please, contact the developer.")
        # end if
        # Get rid of "read_"
        readids_to_seek.append(intern(read_name[5:]))
    # end for

    res_dict = dict()
    i = 0

    for fastq_fpath in fastq_lst:
        fastq_readids = list(get_fastq_readids(fastq_fpath))

        for readid in reversed(readids_to_seek):
            
            if readid in fastq_readids:
                try:
                    res_dict[fastq_fpath] += 1
                except KeyError:
                    res_dict[fastq_fpath] = 1
                finally:
                    readids_to_seek.remove(readid)
                # end try
            # end if
        # end for
        if len(readids_to_seek) == 0:
            break
        # end if
    # end for

    outfile.write("  '{}' -- {} reads;\n".format(fast5_fpath, num_reads))
    for fastq_fpath, curr_nreads in res_dict.items():
        outfile.write("{} reads in '{}';\n".format(curr_nreads, fastq_fpath))
    # end for
    nlost_reads = num_reads - sum(res_dict.values())
    if nlost_reads != 0:
        total_nlost_reads += nlost_reads
        outfile.write("! - {} reads from this FAST5 file are not found in FASTQ files.\n".format(nlost_reads))
    # end if
    outfile.write('-'*20 + '\n')
    print("\r{} - '{}': end processing".format(get_work_time(), os.path.basename(fast5_fpath)))
# end for
outfile.close()

if total_nlost_reads != 0:
    final_verdict = "\n{} reads from FAST5 files are not found in FASTQ files.".format(total_nlost_reads)
# end if
print(final_verdict)

end_time = time()
print( '\n'+get_work_time() + " ({}) ".format(strftime("%Y.%m.%d %H:%M:%S", localtime(end_time))) + "- Task is completed!\n")
platf_depend_exit(0)
