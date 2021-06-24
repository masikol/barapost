#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "3.18.e"
# Year, month, day
__last_update_date__ = "2021-06-24"

# |===== Check python interpreter version =====|

import sys

if sys.version_info.major < 3:
    print( "\nYour python interpreter version is " + "%d.%d" % (sys.version_info.major,
        sys.version_info.minor) )
    print("   Please, use Python 3.\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith("win"):
        raw_input("Press ENTER to exit:")
    # end if
    sys.exit(1)
# end if

from src.platform import platf_depend_exit

# Firstly search for information-providing options (-h and -v):

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print("\n  barapost-local.py\n  Version {}; {} edition;\n".format(__version__, __last_update_date__))
    print("DESCRIPTION:\n")
    print("""This script is designed for taxonomic classification of nucleotide sequences by finding
  the most similar sequence in a nucleotide database stored on local machine.\n""")

    if "--help" in sys.argv[1:]:
        print(""""barapost-local.py" downloads records "discovered" by `barapost-prober.py` (and all replicons
  related to them: other chromosomes, plasmids) from Genbank according to file (`hits_to_download.tsv`)
  generated by "barapost-prober.py" and creates a database on local machine. After that `barapost-local.py classifies
  the rest of data with "BLAST+" toolkit.\n""")
        print("Script processes FASTQ and FASTA (as well as `.fastq.gz` and `.fasta.gz`) files.\n")
        print("""If you have your own FASTA files that can be used as database alone to blast against,
  you can omit "barapost-prober.py" step and go to `barapost-local.py` (see `-l` option).""")
        print("----------------------------------------------------------\n")
        print("Default parameters:\n")
        print("- if no input files are specified, all FASTQ and FASTA files in current directory will be processed;")
        print("- packet size (see `-p` option): 100 sequences;")
        print("- algorithm (see `-a` option): 0 (megaBlast);")
        print("- numbers of threads to launch (`-t` option): 1 thread.")
        print("- BLAST index will be used (see `-i` option).")
        print("----------------------------------------------------------\n")
    # end if

    print("""Files that you want `barapost-local.py` to process should be specified as
  positional arguments (see EXAMPLE #2 running detailed (`--help`) help message).""")
    print("----------------------------------------------------------\n")
    print("OPTIONS:\n")
    print("""-h (--help) --- show help message.
   `-h` -- brief, `--help` -- full;\n""")
    print("-v (--version) --- show version;\n")
    print("""-r (--annot-resdir) --- result directory generated by script `barapost-prober.py`
   This is directory specified to `barapost-prober.py` by `-o` option.
   If you omit `barapost-prober.py` and use your own FASTA files
   to create a database, this directory may not exist before start of `barapost-local.py`
   (i.e. it will be a simple output directory).
   Default value is "barapost_result".\n""")
    print("""-d (--indir) --- directory which contains FASTQ of FASTA files meant to be processed.
   I.e. all FASTQ and FASTA files in this direcory will be processed;\n""")
    print("""-p (--packet-size) --- size of the packet, i.e. number of sequence to align in one blastn launching.
   Value: positive integer number. Default value is 100;\n""")
    print("""-a (--algorithm) --- BLASTn algorithm to use for aligning.
   Available values: 0 for megaBlast, 1 for discoMegablast, 2 for blastn.
   Default is 0 (megaBlast);\n""")
    print("""-l (--local-fasta-to-db) --- your own (local) FASTA file that will be included in
   downloaded database, which "barapost-local.py" creates;\n""")
    print("""-s (--accession) --- accession(s) of GenBank record to download and
   include in database. Multiple accession should be separated by comma without whitespaces;\n""")
    print("-t (--threads) --- number of CPU threads to use;\n")
    print("-i (--use-index) --- whether to use BLAST index to accelerate searches.")
    print("  However, creating index may be memory-consuming.")
    print("  Values: `1` (do use index), `0` (do not use index).\n")

    if "--help" in sys.argv[1:]:
        print("----------------------------------------------------------\n")
        print("EXAMPLES:\n")
        print("""1. Process all FASTA and FASTQ files in working directory with default settings:\n
   barapost-local.py\n""")
        print("""2. Process all files starting with "some_fasta" in the working directory with default settings:\n
   barapost-local.py some_fasta*\n""")
        print("""3. Process one FASTQ file with default settings.
   Directory that contains taxonomic annotation is named `prober_outdir`:\n
   barapost-local.py reads.fastq -r prober_outdir\n""")
        print("""4. Process FASTQ file and FASTA file with discoMegablast, packet size of 50 sequences.
   Directory that contains taxonomic annotation is named `prober_outdir`:\n
   barapost-local.py reads.fastq.gz another_sequences.fasta -a discoMegablast -p 50 -r prober_outdir\n""")
        print("""5. Process all FASTQ and FASTA files in directory named `some_dir`.
   Directory that contains taxonomic annotation is named `prober_outdir`:\n
   barapost-local.py -d some_dir -r prober_outdir\n""")
        print("""6. Process file named `some_reads.fastq`.
   Directory that contains taxonomic annotation is named `prober_outdir`.
   Sequence from file `my_own_sequence.fasta` will be included to the database.
   Launch 4 threads:\n
   barapost-local.py some_reads.fastq -l my_own_sequence.fasta -t 4 -r prober_outdir""")
     # end if
    platf_depend_exit(0)
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if


import os
from glob import glob
from re import search as re_search

import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvd:p:a:r:l:t:s:i:",
        ["help", "version", "indir=", "packet-size=", "algorithm=", "taxannot-resdir=",
        "local-fasta-to-bd=", "threads=", "accession=", "use-index="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

is_fq_or_fa = lambda f: not re_search(r".*\.(m)?f(ast)?(a|q)(\.gz)?$", f) is None

# Default values:
fq_fa_list = list() # list of paths to file meant to be processed
indir_path = None # path to `-d` directory
packet_size = 100
blast_algorithm = "megaBlast"
tax_annot_res_dir = "barapost_result" # directory with taxonomic annotation results
your_own_fasta_lst = list() # list os user's fasta files to be included in database
accs_to_download = list() # list of accessions of GenBank records to download
n_thr = 1 # number of threads
use_index = "true"

# Add positional arguments to fq_fa_list
for arg in args:
    if not is_fq_or_fa(arg):
        print("Error: invalid positional argument: `{}`".format(arg))
        print("Only FAST(A/Q) files can be specified without an option in command line.")
        platf_depend_exit(1)
    # end if
    if not os.path.exists(arg):
        print("Error: file `{}` does not exist!".format(arg))
        platf_depend_exit(1)
    # end if
    fq_fa_list.append( os.path.abspath(arg) )
# end for

for opt, arg in opts:

    if opt in ("-r", "--annot-resdir"):
        tax_annot_res_dir = os.path.abspath(arg)

    elif opt in ("-t", "--threads"):
        try:
            n_thr = int(arg)
            if n_thr < 1:
                raise ValueError
            # end if
        except ValueError:
            print("Error: number of threads must be positive integer number!")
            print(" And here is your value: `{}`".format(arg))
            sys.exit(1)
        # end try
        if n_thr > len(os.sched_getaffinity(0)):
            print("""\nWarning! You have specified {} threads to use
  although {} are available.""".format(n_thr, len(os.sched_getaffinity(0))))
            error = True
            while error:
                reply = input("""\nPress ENTER to switch to {} threads,
  or enter `c` to continue with {} threads,
  or enter `q` to exit:>>""".format(len(os.sched_getaffinity(0)), n_thr))
                reply = reply.lower()
                if reply in ("", 'c', 'q'):
                    error = False
                    if reply == "":
                        n_thr = len(os.sched_getaffinity(0))
                        print("\nNumber of threads switched to {}\n".format(n_thr))
                    elif reply == 'c':
                        pass
                    elif reply == 'q':
                        sys.exit(0)
                    # end if
                else:
                    print("\nInvalid reply!\n")
                # end if
            # end while
        # end if

    elif opt in ("-d", "--indir"):
        if not os.path.isdir(arg):
            print("Error: directory `{}` does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        indir_path = os.path.abspath(arg)

        # Add all fastq and fasta files from `-d` directory to fq_fa_list
        fq_fa_list.extend(list( filter(is_fq_or_fa, glob("{}{}*".format(indir_path, os.sep))) ))

    elif opt in ("-p", "--packet-size"):
        try:
            packet_size = int(arg)
            if packet_size <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("Error: packet_size (`-p` option) must be positive integer number")
            platf_depend_exit(1)
        # end try

    elif opt in ("-l", "--local-fasta-to-bd"):

        if not os.path.exists(arg):
            print("Error: file `{}` does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        your_own_fasta_lst.append(os.path.abspath(arg))

    elif opt in ("-s", "--accession"):

        for acc in arg.split(','):
            accs_to_download.append(acc.partition(".")[0])
        # end for

    elif opt in ("-a", "--algorithm"):
        if not arg in ("0", "1", "2"):
            print("Error: invalid value specified by `-a` option!")
            print("Available values: 0 for megaBlast, 1 for discoMegablast, 2 for blastn")
            print("Your value: `{}`".format(arg))
            platf_depend_exit(1)
        # end if
        blast_algorithm = ("megaBlast", "discoMegablast", "blastn")[int(arg)]

    elif opt in ("-i", "--use-index"):
        if arg in ("0", "1"):
            use_index = "true" if arg == "1" else "false"
        else:
            print("Error: invalid argument passed with `{}` option: `{}`".format(opt, arg))
            print("Available values: `1` (use index) and `0` (do not use index).")
            platf_depend_exit(1)
        # end if
    # end if
# end for


# If no FASTQ or FASTA file have been found
if len(fq_fa_list) == 0:
    # If input directory was specified -- exit
    if not indir_path is None:
        print("""Error: no input FASTQ or FASTA files specified
  or there is no FASTQ and FASTA files in the input directory.\n""")
        platf_depend_exit(1)

    # If input directory was not specified -- look for FASTQ files in working directory
    else:
        fq_fa_list = list(filter( is_fq_or_fa, glob("{}{}*".format(os.getcwd(), os.sep)) ))

        # If there are nothing to process -- just show help message
        if len(fq_fa_list) == 0:
            print("\nbarapost-local.py (Version {})\n".format(__version__))
            print("Usage:")
            print("  barapost-local.py one.fastq.gz another.fasta -r tax_annotation_dir [...] [OPTIONS]")
            print("For more detailed description, run:")
            print("  barapost-local.py -h\n")
            platf_depend_exit(0)
        else:
            # Ask if a user wants to proceed or he/she ran it occasionally and wants just help message
            print("\n {} fasta and/or fastq files are found in working directory.\n".format(len(fq_fa_list)))
            error = True
            while error:
                reply = input("""Press ENTER to process them
  or enter `h` to just see help message:>> """)
                if reply == "":
                    error = False
                elif reply == 'h':
                    error = False
                    print('\n' + '-'*15)
                    print("  barapost-local.py (Version {})\n".format(__version__))
                    print("Usage:")
                    print("  barapost-local.py one.fastq.gz another.fasta -r tax_annotation_dir [...] [OPTIONS]")
                    print("For more detailed description, run:")
                    print("  barapost-local.py -h\n")
                    platf_depend_exit(0)
                else:
                    print("Invalid reply: `{}`\n".format(reply))
                # end if
            # end while
        # end if
    # end if
# end if

# Check if there are duplicated basenames in input files:
for path in fq_fa_list:
    bname = os.path.basename(path)
    same_bnames = tuple(filter(lambda f: os.path.basename(f) == bname, fq_fa_list))
    if len(same_bnames) != 1:
        print("Error: input files must have different base names")
        print("List of files having same name:")
        for p in same_bnames:
            print("`{}`".format(p))
        # end for
        platf_depend_exit(1)
    # end if
# end for
del bname, same_bnames

# Sort list of files that will be processed -- process them in alphabetical order.
fq_fa_list.sort()


# |== Check if 'blast+' tookit is installed ==|

pathdirs = os.environ["PATH"].split(os.pathsep)

# Add '.exe' extention in order to find executables on Windows
if sys.platform.startswith("win"):
    exe_ext = ".exe"
else:
    exe_ext = ""
# end if

for utility in ("blastn"+exe_ext, "makeblastdb"+exe_ext, "makembindex"+exe_ext):

    utility_found = False

    for directory in pathdirs:
        if os.path.exists(directory) and utility in os.listdir(directory):
            utility_found = True
            break
        # end if
    # end for

    if not utility_found:
        print("  Attention!\n`{}` program from BLAST+ toolkit is not installed.".format(utility))
        print("""If this error still occures although you have installed everything
-- make sure that this program is added to PATH variable)""")
        platf_depend_exit(1)
    # end if
# end for

acc_fpath = os.path.join(tax_annot_res_dir, "hits_to_download.tsv")
acc_file_exists = os.path.exists(acc_fpath)

if not acc_file_exists:
    acc_fpath = None
# end if

db_exists = os.path.exists( os.path.join(tax_annot_res_dir, "local_database") )
if db_exists:
    db_exists = db_exists and len(os.listdir(os.path.join(tax_annot_res_dir, "local_database"))) != 0
# end if

# Create result directory, if it does not exist
if not os.path.exists(tax_annot_res_dir):
    try:
        os.makedirs(tax_annot_res_dir)
    except OSError as oserr:
        print("Error: unable to create output directory.")
        print( str(oserr) )
        print("Prober just tried to create directory `{}` and crashed.".format(tax_annot_res_dir))
        platf_depend_exit(1)
    # end try
# end if

if not db_exists:

    # Check if there is a way to create a database (or to use old one):
    if not acc_file_exists and len(your_own_fasta_lst) == 0 and len(accs_to_download) == 0:

        print("Error: cannot create a database. Reasons:")
        print("1. There is no `hits_to_download.tsv` file in directory `{}`.".format(tax_annot_res_dir,))
        print("2. No fasta file is passed with `-l` option and not accession number(s) are passed with `-s` option")
        if tax_annot_res_dir == "barapost_result":
            print("""\nMaybe, the reason is that output directory generated by `barapost-prober.py`
      isn't named `barapost_result` and you have forgotten to specify it with `-r` option?""")
        # end if
        platf_depend_exit(1)
    # end if
# end if
del db_exists, acc_file_exists

# Generate path to log file:
from src.platform import get_logfile_path

from src.printlog import get_full_time, printn, printlog_info, printlog_info_time
from src.printlog import printlog_error, printlog_error_time, log_info, printlog_warning
import logging
logging.basicConfig(filename=get_logfile_path("barapost-local", tax_annot_res_dir),
    format='%(levelname)s: %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO, filemode='w')
log_info(sys.platform)
log_info(sys.implementation)
log_info(sys.version)

print("|=== barapost-local.py (version {}) ===|\n".format(__version__))
log_info("barapost-local.py (version {})".format(__version__))
print(get_full_time() + "- Start working\n")
log_info("Start working.")


#                       |===== Proceed =====|

printlog_info(" - Logging to `{}`".format(logging.getLoggerClass().root.handlers[0].baseFilename))
printlog_info(" - Output directory: `{}`;".format(tax_annot_res_dir))
printlog_info(" - Packet size: {} sequences;".format(packet_size))
printlog_info(" - BLAST algorithm: {};".format(blast_algorithm))
printlog_info(" - Threads: {};".format(n_thr))
print()

s_letter = '' if len(fq_fa_list) == 1 else 's'
printlog_info(" {} file{} will be processed.".format( len(fq_fa_list), s_letter))
if len(fq_fa_list) != 1:
    log_info("Here they are:")
else:
    log_info("Here it is:")
# end if
for i, path in enumerate(fq_fa_list):
    log_info("  {}. `{}`".format(i+1, path))
# end for

if len(accs_to_download) != 0:
    print()
    printlog_info("Following GenBank records were specified in command line \
and will be downloaded and included in database:")
    for i, acc in enumerate(accs_to_download):
        printlog_info("  {}. {}".format(i+1, acc))
# end if

if len(your_own_fasta_lst) != 0:
    preposition = " besides downloaded ones" if not acc_fpath is None else ""
    print()
    printlog_info("Following FASTA files will be included in database{}:".format(preposition))
    for i, path in enumerate(your_own_fasta_lst):
        printlog_info("  {}. `{}`".format(i+1, path))
    # end for
# end if

printlog_info('-'*30)

# Algorithms in 'blast+' are named in a little different way comparing to BLAST server.
# In order to provide full interface compatibility with 'barapost-prober.py' we will merely change values here.
if blast_algorithm == "megaBlast":
    blast_algorithm = "megablast"
elif blast_algorithm == "discoMegablast":
    blast_algorithm = "dc-megablast"
# end if

import src.legacy_taxonomy_handling as legacy_taxonomy_handling

# Form path to taxonomy file:
taxonomy_dir = os.path.join(tax_annot_res_dir, "taxonomy")
if not os.path.isdir(taxonomy_dir):
    try:
        os.makedirs(taxonomy_dir)
    except OSError as err:
        printlog_error_time("Error: cannot create taxonomy directory `{}`".format(taxonomy_dir))
        printlog_error_time(str(err))
        platf_depend_exit(1)
    # end try
# end if
taxonomy_path = os.path.join(taxonomy_dir, "taxonomy.tsv")

# Check if there is legacy taxonomy file and, if so, reformat it to new (TSV) format
legacy_taxonomy_handling.check_deprecated_taxonomy(tax_annot_res_dir)

from src.barapost_local_modules.build_local_db import build_local_db

# Indexed discontiguous searches are not supported:
#    https://www.ncbi.nlm.nih.gov/books/NBK279668/#usermanual.Megablast_indexed_searches
if use_index == "true" and blast_algorithm == "dc-megablast":
    printlog_warning("Warning: BLAST index cannot be used in alignment algorithm is DiscoMegablast.")
    printlog_warning("Index will be created anyway.")
# end if

# Build a database
db_path = build_local_db(tax_annot_res_dir,
    acc_fpath,
    your_own_fasta_lst,
    accs_to_download,
    use_index)

if blast_algorithm == "dc-megablast":
    use_index = "false"
# end if

if use_index == "true" and len(glob(os.path.join(tax_annot_res_dir, "local_database", "*idx"))) == 0:
    printlog_warning("Warning: unable to use BLAST index: `*.idx` files don't exist in output directory.")
    printlog_warning("Try to rebuild the database next time.")
    printlog_warning("Now, BLAST searches just will not use index.")
    use_index = "false"
# end if


# Create temporary directory for query files:
queries_tmp_dir = os.path.join(tax_annot_res_dir, "queries-tmp")
if not os.path.isdir(queries_tmp_dir):
    try:
        os.makedirs(queries_tmp_dir)
    except OSError as oserr:
        printlog_error_time("unable to create directory for queries:")
        printlog_error("  `{}`".format(queries_tmp_dir))
        printlog_error("Reason: {}".format(str(oserr) ))
        platf_depend_exit(1)
    # end try
# end if


# Proceeding.
# The main goal of multiprocessing is to isolate processes from one another.
#
# Two situations are available:
#   1. Number of threads <= number of files meant to be processed ('many_files'-parallel mode):
#      Files will be distribured equally among processes.
#      Processes interact with one another only while printing things to the console
#      for user's entertainment.
#   2. Number of threads > number of files meant to be processed ('few_files'-parallel mode):
#      Files will be processed one by one. They will be divided into equal blocks,
#      and these blocks will be distributed among processes.
#      Processes interact with one another while writing to result file and
#      while printing things to the console.

print()
printlog_info_time("Starting classification.")
printn("  Working...")

if n_thr <= len(fq_fa_list):
    if n_thr != 1:

        # Proceed 'many_files'-parallel processing

        from src.barapost_local_modules.parallel_mult_files import process

        process(fq_fa_list,
            n_thr,
            packet_size,
            tax_annot_res_dir,
            blast_algorithm,
            use_index,
            db_path)

    else:

        # Proceed single-thread processing

        from src.barapost_local_modules.single_thread_mult_files import process

        process(fq_fa_list,
            packet_size,
            tax_annot_res_dir,
            blast_algorithm,
            use_index,
            db_path)
    # end if
else:

    # Proceed 'few_files'-parallel processing

    from src.barapost_local_modules.parallel_single_file import process

    process(fq_fa_list,
        n_thr,
        packet_size,
        tax_annot_res_dir,
        blast_algorithm,
        use_index,
        db_path)
# end if

# Remove everything in 'queries_tmp_dir'
try:
    for qpath in glob( os.path.join(queries_tmp_dir, '*') ):
        os.unlink(qpath)
    # end for
    os.rmdir(queries_tmp_dir)
except OSError as oserr:
    printlog_error_time("unable to remove directory `{}`".format(queries_tmp_dir))
    printlog_error( str(oserr) )
    printlog_error("Don't worry -- barapost-local has completed it's job just fine,")
    printlog_error("   the only thing that some temporary files are left in the directory mentioned above.\n")
# end try

print("\r{}".format(' ' * len("Working...")))
print("{} - Task is completed!".format(get_full_time()))
log_info("Task is completed!")
platf_depend_exit(0)
