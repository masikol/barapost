#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.1.d"
# Year, month, day
__last_update_date__ = "2020-12-28"

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
            with open(outfpath, 'w') as outfile:
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

print( get_work_time() + " ({}) ".format(strftime("%Y-%m-%d %H:%M:%S", localtime(start_time))) + "- Start working\n")

print("{} FAST5 files will be processed.".format(len(fast5_lst)))
print("""Results will be written to the following file:
  '{}'\n""".format(os.path.abspath(outfpath)))

print("{} - Primary validation...".format(get_work_time()))

from threading import Thread

# Value that contains number of processed files:
global inc_val
# Flag that signals printer to stop
global stop

def status_printer(get_inc_val):
    """
    Function meant to be launched as threading.Thread in order to indicate progress each second.

    :param get_inc_val: function that returns 'inc_val' value -- the number of processed files;
    """
    printn("{} - 0/{} files processed. Working...".format(get_work_time(), len(fast5_lst)))
    saved_val = get_inc_val()
    while not stop:
        loc_inc_val = get_inc_val()
        printn("\r{} - {}/{} files processed. Working...".format(get_work_time(), loc_inc_val, len(fast5_lst)))
        if loc_inc_val != saved_val:
            saved_val = get_inc_val()
        sleep(1)
    # end while
    printn("\r{} - {}/{} files processed.".format(get_work_time(), get_inc_val(), len(fast5_lst)) +
        ' '*len(" Working..."))
# end def status_printer

inc_val = 0
get_inc_val = lambda: inc_val

printer = Thread(target=status_printer, args=(get_inc_val,), daemon=True) # create thread
stop = False # raise the flag
printer.start() # start waiting

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
        inc_val += 1
    # end try
# end for

stop = True # lower the flag
printer.join()
print("\n{} - Primary validation -- Ok\n".format(get_work_time()))

from gzip import open as open_as_gzip
is_gzipped = lambda f: True if f.endswith(".gz") else False

OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line,   # format text line
    lambda line: line.decode("utf-8")  # format gzipped line
)


# According to
# https://github.com/nanoporetech/ont_h5_validator/blob/master/h5_validator/schemas/multi_read_fast5.yaml
ont_read_signature = r"([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})"


def fmt_read_id(read_id):

    srch_ont_read = re_search(ont_read_signature, read_id)
    if srch_ont_read is None:
        return read_id.partition(' ')[0]
    else:
        return srch_ont_read.group(1)
# end def fmt_read_id


def get_fastq_readids(fastq_fpath):

    how_to_open = OPEN_FUNCS[ is_gzipped(fastq_fpath) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq_fpath) ]

    with how_to_open(fastq_fpath) as fastq_file:
        line = fmt_func(fastq_file.readline())
        yield fmt_read_id(line)
        i = 1
        while line != "":
            line = fmt_func(fastq_file.readline())
            i += 1
            if i == 4:
                i = 0
            elif i == 1:
                yield fmt_read_id(line)
            # end if
        # end while
    # end with
# end def get_fastq_readids

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


final_verdict = "\nAll reads from all FAST5 files are found"
total_nlost_reads = 0

inc_val = 0
get_inc_val = lambda: inc_val
print("Searching for reads from FAST5 files started")

printer = Thread(target=status_printer, args=(get_inc_val,), daemon=True) # create thread
stop = False # raise the flag
printer.start() # start waiting

outfile = open(outfpath, 'w')
for fast5_fpath in fast5_lst:

    f5_file = h5py.File(fast5_fpath, 'r')

    readids_to_seek = fast5_readids(f5_file) # list of not-found-yet read IDs
    readids_to_seek = list(map(fmt_read_id, readids_to_seek))
    num_reads = len(readids_to_seek)
    res_dict = dict()

    for fastq_fpath in fastq_lst:
        fastq_readids = list(get_fastq_readids(fastq_fpath))

        for readid in readids_to_seek:
            
            if readid in fastq_readids:
                try:
                    res_dict[fastq_fpath] += 1
                except KeyError:
                    res_dict[fastq_fpath] = 1
                # end try
            # end if
        # end for
        if len(readids_to_seek) == 0:
            break
        # end if
    # end for

    inc_val += 1

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
    outfile.flush()

# end for
outfile.close()

stop = True # lower the flag
printer.join()

if total_nlost_reads != 0:
    final_verdict = "\n{} reads from FAST5 files are not found in FASTQ files.".format(total_nlost_reads)
# end if
print(final_verdict)

end_time = time()
print( '\n'+get_work_time() + " ({}) ".format(strftime("%Y-%m-%d %H:%M:%S", localtime(end_time))) + " - Task is completed!\n")
platf_depend_exit(0)
