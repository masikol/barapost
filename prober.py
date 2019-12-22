#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.14.a"
# Year, month, day
__last_update_date__ = "2019-12-22"

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
    Function asks to press ENTER press on Windows
        and exits after that.

    :type exit_code: int;
    """
    if platform.startswith("win"):
        input("Press ENTER to exit:")
    # end if
    exit(exit_code)
# end def platf_depend_exit

from sys import argv

# First search for information-providing options:

if "-h" in argv[1:] or "--help" in argv[1:]:

    print("\n  prober.py\n  Version {}; {} edition;\n".format(__version__, __last_update_date__))
    print("DESCRIPTION:\n")
    print("""This script is designed for determinating the taxonomic position
  of nucleotide sequences by sending each of them to NCBI BLAST server and regarding the best hit.\n""")

    if "--help" in argv[1:]:
        print("""The main goal of this script is to send a probing batch (see `-b` option) of sequences to
      NCBI BLAST server and discover, what Genbank records can be downloaded and used for
      building a database on your local machine by "barapost.py".\n""")
        print("""This script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.\n
    Results of taxonomic annotation are written in TSV file named according to name of input file(s),
      that can be found in result directory.\n""")
        print(""""prober.py" also generates a file named 'hits_to_download.tsv'. It contains accessions and
      names of Genbank records that can be used for building a database on your local machine by "barapost.py".""")
        print("----------------------------------------------------------\n")
        print("Default parameters:\n")
        print(" - all FASTQ and FASTA files in current directory will be processed;")
        print(" - packet size (see '-p' option): 100 sequences;")
        print(" - probing batch size (see '-b' option): 200 sequences;")
        print(" - algorithm (see '-a' option): `megaBlast`;")
        print(" - organisms (see '-g' option): full 'nt' database, i.e. no slices;")
        print(" - output directory ('-o' option): directory named 'barapost_result'")
        print("   nested in current directory;")
        print(" - no email information ('-e' option) is send to NCBI;")
        print(" - prober sends sequences intact (i.e. does not prune them (see '-x' option));")
        print("----------------------------------------------------------\n")
    # end if

    print("""Files that you want 'prober.py' to process should be specified as
  positional arguments (see EXAMPLE #2 running detailed (--help) help message).
  Wildcards do work: './prober.py my_directory/*' will process all files in 'my_directory'.""")
    print("----------------------------------------------------------\n")
    print("OPTIONS:\n")
    print("""-h (--help) --- show help message.
        '-h' -- brief, '--help' -- full;\n""")
    print("-v (--version) --- show version;\n")
    print("""-d (--indir) --- directory which contains FASTQ of FASTA files meant to be processed.
        I.e. all FASTQ and FASTA files in this direcory will be processed;
        Input files can be gzipped.\n""")
    print("""-o (--outdir) --- output directory.
        Default value: 'barapost_result';\n""")
    print("""-p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
        Value: integer number [1, 500]. Default value is 100;\n""")
    print("""-a (--algorithm) --- BLASTn algorithm to use for aligning.
        Available values: 'megaBlast', 'discoMegablast', 'blastn'.
        Default is megaBlast;\n""")
    print("""-g (--organisms) --- TaxIDs of organisms to align your sequences against. I.e. 'nt' database slices.
        More clearly, functionality of this option is totally equal to "Organism" text boxes
        on this BLASTn page:
         'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'.
        Format of value: 
          <organism1_taxid>,<organism2_taxid>...
        See EXAMPLES #2 and #3 below.
        Spaces are not allowed.
        Default value is full 'nt' database, i.e. no slices.
        You can find your Taxonomy IDs here: 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi'.\n""")
    print("""-b (--probing-batch-size) --- total number of sequences that will be sent to BLAST server
        during 'prober.py' run.
        You can specify '-b all' to process all your sequeces by 'prober.py'.
        Value: positive integer number.
        Default value is 200;\n""")
    print("""-e (--email) --- your email. Please, specify your email when you run "prober.py",
        so that the NCBI can contact you if there is a problem. See EXAMPLE #2 below.\n""")
    print("""-x (--max-seq-len) --- maximum length of a sequence that will be sent to NCBI BLAST.
        It means that prober can prune your sequences before sending in order to spare NCBI servers.
        This feature is disabled by default;""")

    if "--help" in argv[1:]:
        print("----------------------------------------------------------\n")
        print("EXAMPLES:\n")
        print("""1. Process all FASTA and FASTQ files in working directory with default settings:\n
  ./prober.py\n""")
        print("""2. Process all files in the working directory that start with "some_my_fasta".
  Provide NCBI with your email. Use default settings:\n
  ./prober.py some_my_fasta* -e my.email@smth.com\n""")
        print("""3. Process one file with default settings:\n
  ./prober.py reads.fastq\n""")
        print("""4. Process a FASTQ file and a FASTA file with discoMegablast, packet size of 100 sequences.
        Search only among Erwinia sequences (551 is Erwinia taxid):\n
  ./prober.py reads_1.fastq.gz some_sequences.fasta -a discoMegablast -p 100 -g 551\n""")
        print("""5. Process all FASTQ and FASTA files in directory named `some_dir`.
  Process 300 sequences, packet size is 100 sequnces (3 packets will be sent).
  Search only among Escherichia (taxid 561) and viral (taxid 10239) sequences:\n
  ./prober.py -d some_dir -g 561,10239 -o outdir -b 300 -p 100""")
    # end if
    platf_depend_exit(0)
# end if

if "-v" in argv[1:] or "--version" in argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, sleep, gmtime
start_time = time()

def get_work_time():
    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def get_work_time

# |===========================================|

import os
from re import search as re_search
from re import match as re_match
from glob import glob

def err_fmt(text):
    """Function for configuring error messages"""
    return "\n   \a!! - ERROR: " + text + '\n'
# end def print_error


from sys import stdout as sys_stdout
def printn(text):
    """
    Function prints text to the console without adding '\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
    """
    sys_stdout.write(text)
    sys_stdout.flush() # display text immediately
# end def printn

import getopt

try:
    opts, args = getopt.gnu_getopt(argv[1:], "hvd:o:p:a:g:b:e:x:",
        ["help", "version", "indir=", "outdir=", "packet-size=",
        "algorithm=", "organisms=", "probing-batch-size=", "email=",
        "max-seq-len="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

is_fq_or_fa = lambda f: True if not re_search(r".*\.(m)?f(ast)?(a|q)(\.gz)?$", f) is None else False

# Default values:
fq_fa_list = list()
indir_path = None
outdir_path = "barapost_result"
packet_size = 100
probing_batch_size = 200
send_all = False
blast_algorithm = "megaBlast"
taxid_list = list() # default is whole 'nt' database
user_email = ""
max_seq_len = None # maximum length of a sequence sent to NCBI

# Add positional arguments to fq_fa_list
for arg in args:
    if not os.path.exists(arg) or not is_fq_or_fa(arg):
        print(err_fmt("invalid positional argument: '{}'".format(arg)))
        print("Only FAST(A/Q) files can be specified without an option in command line.")
        platf_depend_exit(1)
    # end if
    fq_fa_list.append( os.path.abspath(arg) )
# end for

for opt, arg in opts:

    if opt in ("-d", "--indir"):
        if not os.path.exists(arg):
            print(err_fmt("directory '{}' does not exist!".format(arg)))
            platf_depend_exit(1)
        # end if
        
        if not os.path.isdir(arg):
            print(err_fmt("'{}' is not a directory!".format(arg)))
            platf_depend_exit(1)
        # end if
        
        indir_path = os.path.abspath(arg)

        fq_fa_list.extend(list( filter(is_fq_or_fa, glob("{}{}*".format(indir_path, os.sep))) ))

    elif opt in ("-o", "--outdir"):
        outdir_path = os.path.abspath(arg)

    elif opt in ("-p", "--packet-size"):
        try:
            packet_size = int(arg)
            if packet_size < 1 or packet_size > 500:
                raise ValueError
            # end if
        except ValueError:
            print(err_fmt("packet_size (-p option) must be integer number from 1 to 500"))
            platf_depend_exit(1)
        # end try

    elif opt in ("-a", "--algorithm"):
        if not arg in ("megaBlast", "discoMegablast", "blastn"):
            print(err_fmt("invalid value specified by '-a' option!"))
            print("Available values: 'megaBlast', 'discoMegablast', 'blastn'")
            platf_depend_exit(1)
        # end if
        blast_algorithm = arg

    elif opt in ("-g", "--organisms"):

        taxid_list = arg.strip().split(',')

        try:
            for taxid in taxid_list:
                buff_var = int(taxid)
                if buff_var < 0:
                    raise ValueError
                # end if
            # end for
        except ValueError:
            print(err_fmt("TaxID should be positive integer number\a"))
            platf_depend_exit(1)
        # end try

    elif opt in ("-b", "--probing-batch-size"):
        # Switch 'send_all' to True in order to process all sequences
        if arg == "all":
            send_all = True
            continue
        # end if
        try:
            probing_batch_size = int(arg)
            if probing_batch_size <= 0:
                raise ValueError
            # end if
        except ValueError:
            print(err_fmt("probing batch size ('-b' option) must be positive integer number!"))
            platf_depend_exit(1)
        # end try

    elif opt in ("-e", "--email"):
        if arg != "" and re_match(r".+@.+\..+", arg) is None:
            print("Your email does not seem like an email.")
            print("Please check it and, if you are sure that your email is right, \
\n   but prober still refuses it  -- contact the developer.")
            platf_depend_exit(1)
        # end if
        user_email = arg

    elif opt in ("-x", "--max-seq-len"):

        try:
            max_seq_len = int(arg)
            if max_seq_len <= 0:
                raise ValueError
            # end if
        except ValueError:
            print(err_fmt("maximum sequence length must bu positive interger number!"))
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try
    # end if
# end for


# If no FASTQ or FASTA file have been found
if len(fq_fa_list) == 0:
    # If input directory was specified -- exit
    if not indir_path is None:
        print(err_fmt("""no input FASTQ or FASTA files specified
    or there is no FASTQ and FASTA files in the input directory.\n"""))
        platf_depend_exit(1)
    
    # If input directory was not specified -- look for FASTQ files in current directory
    else:
        fq_fa_list = list(filter( is_fq_or_fa, glob("{}{}*".format(os.getcwd(), os.sep)) ))
        if len(fq_fa_list) == 0:
            print(err_fmt("there is no FASTQ or FASTA files to process found."))
            platf_depend_exit(1)
        # end if
    # end if
# end if

# Sort list of files that will be processed -- process them in alphabetical order.
fq_fa_list.sort()

from gzip import open as open_as_gzip # input files might be gzipped
from xml.etree import ElementTree # for retrieving information from XML BLAST report
from sys import intern

import http.client
import urllib.request
from urllib.error import HTTPError
import urllib.parse
import socket
import shelve

# |===== Functionality for proper processing of gzipped files =====|

OPEN_FUNCS = (open, open_as_gzip)

is_gzipped = lambda file: True if file.endswith(".gz") else False
is_fastq = lambda f: True if not re_search(r".*\.fastq(\.gz)?$", f) is None else False

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)

# Delimiter for result tsv file:
DELIM = '\t'

# File format constants:
FASTQ_LINES_PER_READ = 4
FASTA_LINES_PER_SEQ = 2

taxonomy_dir = os.path.join(outdir_path, "taxonomy")
if not os.path.isdir(taxonomy_dir):
    os.makedirs(taxonomy_dir)
# end if
indsxml_path = os.path.join(taxonomy_dir, "indsxml.gbc.xml")
taxonomy_path = os.path.join(taxonomy_dir, "taxonomy")


# |=== Check if there are enough sequeneces in files (>= probing_batch_size) ===|
seqs_at_all = 0
print() # just print new line
for file in fq_fa_list:
    how_to_open = OPEN_FUNCS[ is_gzipped(file) ]
    if is_fastq(file):
        seqs_at_all += int( sum(1 for line in how_to_open(file)) / FASTQ_LINES_PER_READ )
    else:
        fmt_func = FORMATTING_FUNCS[ is_gzipped(file) ]

        seqs_at_all = len(tuple(filter(lambda l: True if l.startswith('>') else False,
            map(fmt_func, how_to_open(file).readlines()))))
    # end if
    if seqs_at_all >= probing_batch_size:
        break
    # end if
# end for


# Print a warning message if a user has specified batch size that is greater than number of sequences he has at all.
# And do not disturb him if he has run 'prober.py' with default batch size.
if seqs_at_all < probing_batch_size and ("-b" in argv or "--probing_batch_size" in argv):
    if send_all:
        probing_batch_size = seqs_at_all
    else:
        print('\n'+'-'*20)
        print("\a  Warning!\n There are totally {} sequences in your files.".format(seqs_at_all))
        print(" Probing batch size specified by you is {}".format(probing_batch_size))

        while True:
            reply = input("\nPress ENTER to process all your sequences anyway.\n  Or enter 'q' to exit:>>")
            if reply == "":
                probing_batch_size = seqs_at_all
                break
            elif reply == 'q':
                platf_depend_exit(0)
            else:
                print(err_fmt("invalid reply"))
                continue
            # end if
        # end while
    # end if
# end if

if seqs_at_all < probing_batch_size and not ("-b" in argv or "--probing_batch_size" in argv):
    probing_batch_size = seqs_at_all
# end if

packet_size = min(packet_size, probing_batch_size)

if not os.path.isdir(outdir_path):
    try:
        os.makedirs(outdir_path)
    except OSError as oserr:
        print_error("unable to create result directory")
        print( str(oserr) )
        print("Prober just tried to create directory '{}' and crushed.".format(outdir_path))
        platf_depend_exit(1)
    # end try
# end if

# |===== Function for checking if 'https://blast.ncbi.nlm.nih.gov' is available =====|

def check_connection():
    """
    Function checks if 'https://blast.ncbi.nlm.nih.gov' is available.

    :return: None if 'https://blast.ncbi.nlm.nih.gov' is available;
    """
    printn("Checking Internet connection...")

    if not platform.startswith("win"):
        check_mark = "\u2714"
    else:
        check_mark = "OK"
    # end if

    try:
        ncbi_server = "https://blast.ncbi.nlm.nih.gov"
        status_code = urllib.request.urlopen(ncbi_server).getcode()
        # Just in case
        if status_code != 200:
            print('\n' + get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            print("Check your Internet connection.\a")
            print("Status code: {}".format(status_code))
            platf_depend_exit(-2)
        # end if
    except OSError as err:
        print('\n' + get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
        print("Check your Internet connection.\a")
        print( str(err) )

        # 'urllib.request.HTTPError' can provide a user with information about the error
        if isinstance(err, HTTPError):
            print("Status code: {}".format(err.code))
            print(err.reason)
        # end if
        platf_depend_exit(-2)
    else:
        print("\rChecking Internet connection... {}".format(check_mark))
    # end try
# end def check_connection

check_connection()

# There some troubles with file extention on Windows, so let's make a .txt file for it:
log_ext = ".log" if not platform.startswith("win") else ".txt"
logfile_path = os.path.join(outdir_path, "prober_log_{}{}".format(strftime("%Y-%m-%d_%H-%M-%S", localtime(start_time)), log_ext))
logfile = open(logfile_path, 'w')

def printl(text=""):
    """
    Function for printing text to console and to log file.
    """
    print(text)
    logfile.write(str(text).strip('\r') + '\n')
    logfile.flush()
# end def printl

def println(text=""):
    """
    Function for printing text to console and to log file.
    The only difference from 'printl' -- text that is printed to console does not end with '\\n'
    """
    printn(text)
    logfile.write(str(text).strip('\r') + '\n')
    logfile.flush()
# end def printl

printl("\n |=== prober.py (version {}) ===|\n".format(__version__))
printl( get_work_time() + " ({}) ".format(strftime("%Y-%m-%d %H:%M:%S", localtime(start_time))) + "- Start working\n")


def verify_taxids(taxid_list):
    """
    Funciton verifies TaxIDs passed to 'prober' with -g option.
    Function requests NCBI Taxonomy Browser and parses organism's name from HTML response.
    What is more, this function configures 'oraganisms' list - it will be included to BLAST requests.

    :param taxid_list: list of TaxIDs. TaxIDs are strings, but they are verified to be integers
        during CL argument parsing;
    :type taxid_list: list<str>;

    Returns list of strings of the following format: "<tax_name> (taxid:<TaxID>)>"
    """
    organisms = list()
    if len(taxid_list) > 0:

        printl("Verifying TaxIDs:")
        if not platform.startswith("win"):
            check_mark = " \u2714"
            cross = " \u274c"
        else:
            check_mark = ""
            cross = ""
        # end if
        for taxid in taxid_list:
            println("   {} - ".format(taxid))
            tax_url = "https://ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}".format(taxid)
            try:
                tax_resp = urllib.request.urlopen(tax_url)
                tax_name = re_search(r"Taxonomy browser \((.+?)\)", tax_resp.read().decode("utf-8")).group(1)
            except AttributeError:
                println("\aError: TaxID not found")
                print(cross)
                print("Please, check your TaxID: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi")
                platf_depend_exit(1)
            except OSError as oserr:
                printl(err_fmt("something is wrong with connection:"))
                printl( str(oserr) )
                platf_depend_exit(-2)
            else:
                println(tax_name)
                print(check_mark)
                organisms.append("{} (taxid:{})".format(tax_name, taxid))
            # end try
        # end for
        printl('-'*30 + '\n')

    # end if
    return organisms
# end def verify taxids


def ask_for_resumption():
    """
    Function asks a user if he/she wants to resume the previous run.

    Returns True if the decision is to resume, else False;
    """

    resume = None

    while resume is None:
        resume = input("""
Would you like to resume the previous run?
    1. Resume!
    2. Start from the beginning.

Enter the number (1 or 2):>> """)
        # Check if entered value is integer number. If no, give another attempt.
        try:
            resume = int(resume)
            # Check if input number is 1 or 2
            if resume != 1 and resume != 2:
                print("\n   Not a VALID number entered!\a\n" + '~'*20)
                resume = None
            else:
                print("You have chosen number " + str(resume) + '\n')
            # end if
        except ValueError:
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            resume = None
        # end try

    return True if resume == 1 else False
# end def ask_for_resumption


# |===== End of question funtions =====|

get_phred33 = lambda q_symb: ord(q_symb) - 33

def get_read_avg_qual(qual_str):

    phred33 = map(get_phred33, list(qual_str))
    read_qual = round( sum(phred33) / len(qual_str), 2 )
    return read_qual
# end def get_read_avg_qual


def configure_qual_dict(fastq_path):

    qual_dict = dict()
    how_to_open = OPEN_FUNCS[ is_gzipped(fastq_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq_path) ]

    with how_to_open(fastq_path) as fastq_file:
        counter = 1
        line = fmt_func(fastq_file.readline())
        while line != "":
            if counter == 1:
                seq_id = intern( fmt_read_id(line).replace('@', '') )
            # end if
            
            counter += 1
            line = fmt_func(fastq_file.readline())
            if counter == 4:
                qual_dict[seq_id] = get_read_avg_qual(line)
                counter = 0
            # end if
        # end while
    # end with

    return qual_dict
# end def configure_qual_dict


from shutil import copyfileobj as shutil_copyfileobj

def fastq2fasta(fq_fa_path, i, new_dpath):
    """
    Function conwerts FASTQ file to FASTA format, if there is no FASTA file with
        the same name as FASTQ file. Also it counts sequences in this file.

    :param fq_fa_path: path to FASTQ or FASTA file being processed;
    :type fq_fa_path: str;
    :param i: order number of fq_fa_path;
    :type i: int;
    :param new_dpath: path to current (corresponding to fq_fa_path file) result directory;
    :type new_dpath: str;

    Returns dict of the following structure:
    {
        "fpath": path_to_FASTA_file (str),
        "nreads": number_of_reads_in_this_FASTA_file (int)
    }
    """
    
    fasta_path = re_search(r"(.*)\.(m)?f(ast)?(a|q)", os.path.basename(fq_fa_path)).group(1) + ".fasta"
    fasta_path = os.path.join(new_dpath, fasta_path) # place FASTA file into result directory

    how_to_open = OPEN_FUNCS[ is_gzipped(fq_fa_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_fa_path) ]

    fastq_patt = r".*\.f(ast)?q(\.gz)?$"

    num_lines = 0 # variable for counting lines in a file
    if not re_search(fastq_patt, fq_fa_path) is None and not os.path.exists(fasta_path+".gz"):

        printl("\n{}. '{}' --> FASTA".format(i+1, os.path.basename(fq_fa_path)))

        global FASTQ_LINES_PER_READ

        with how_to_open(fq_fa_path) as fastq_file, open(fasta_path, 'w') as fasta_file:

            counter = 1 # variable for retrieving only 1-st and 2-nd line of FASTQ record
            for line in fastq_file:
                line = fmt_func(line)
                # write only 1-st and 2-nd line out of 4
                if counter <= 2:
                    if line[0] == '@':
                        line = '>' + line[1:]  # replace '@' with '>'
                    # end if
                    fasta_file.write(line + '\n')
                # reset the counter if the 4-th (quality) line has been encountered
                elif counter == 4:
                    counter = 0
                # end if
                counter += 1
                num_lines += 1
            # end for
        # end with
        num_reads = int(num_lines / FASTQ_LINES_PER_READ) # get number of sequences

        print("Gzipping '{}'...".format(fasta_path))

        # GNU gzip utility is faster, but there can be presence of absence of it
        gzip_util = "gzip"
        util_found = False
        for directory in os.environ["PATH"].split(os.pathsep):
            if os.path.isdir(directory) and gzip_util in os.listdir(directory):
                util_found = True
                break
            # end if
        # end for

        if util_found:
            os.system("{} {}".format(gzip_util, fasta_path))
        else:
            # form .fasta.gz file 'by hand'
            with open(fasta_path, 'rb') as fasta_file, open_as_gzip(fasta_path+".gz", "wb") as fagz_file:
                shutil_copyfileobj(fasta_file, fagz_file)
            # end with
            os.unlink(fasta_path) # remove plain FASTA file
        # end if
    
    # IF FASTA file is already created
    # We need only number of sequences in it.
    elif not re_search(fastq_patt, fq_fa_path) is None and os.path.exists(fasta_path+".gz"):
        num_lines = sum(1 for line in how_to_open(fq_fa_path)) # get number of lines
        num_reads = int( num_lines / FASTQ_LINES_PER_READ ) # get number of sequences
    # We've got FASTA source file
    # We need only number of sequences in it.
    else:

        num_reads = len(tuple(filter(lambda l: True if l.startswith('>') else False,
            map(fmt_func, how_to_open(fq_fa_path).readlines()))))
        fasta_path = fq_fa_path
        return {"fpath": fasta_path, "nreads": num_reads}
    # end if

    return {"fpath": fasta_path+".gz", "nreads": num_reads}
# end def fastq2fasta



def rename_file_verbosely(file, directory):

    if os.path.exists(file):

        is_analog = lambda f: file[:os.path.basename(file).rfind('.')] in f
        num_analog_files = len( list(filter(is_analog, os.listdir(directory))) )

        printl('\n' + get_work_time() + " - Renaming old file:")
        name_itself = file[: file.rfind('.')]
        ext = file[file.rfind('.'):]
        num_analog_files = str(num_analog_files)
        new_name = name_itself+"_old_"+num_analog_files+ext

        printl("   '{}' --> '{}'".format(os.path.basename(file), new_name))
        os.rename(file, new_name)
    # end if
# end def rename_file_verbosely


def look_around(outdir_path, new_dpath, fasta_path, blast_algorithm):
    """
    Function looks around in order to ckeck if there are results from previous runs of this script.

    Returns None if there is no result from previous run.
    If there are results from previous run, returns a dict of the following structure:
    {
        "pack_size": packet_size (int),
        "sv_npck": saved_number_of_sent_packet (int),
        "RID": saved_RID (str),
        "tsv_respath": path_to_tsv_file_from_previous_run (str),
        "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file (int),
        "tmp_fpath": path_to_pemporary_file (str)
    }

    :param new_dpath: path to current (corresponding to fq_fa_path file) result directory;
    :type new_dpath: str;
    :param fasta_path: path to current (corresponding to fq_fa_path file) FASTA file;
    :type fasta_path: str;
    :param blast_algorithm: BLASTn algorithm to use.
        This parameter is necessary because it is included in name of result files;
    :type blast_algorithm: str;
    """

    # "hname" means human readable name (i.e. without file path and extention)
    fasta_hname = os.path.basename(fasta_path) # get rid of absolute path
    fasta_hname = re_search(r"(.*)\.(m)?f(ast)?a", fasta_hname).group(1) # get rid of '.fasta' extention
    how_to_open = OPEN_FUNCS[ is_gzipped(fasta_path) ]

    # Form path to temporary file
    tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)
    # Form path to result file
    tsv_res_fpath = "{}.tsv".format(os.path.join(new_dpath, "classification"))
    # Form path to accession file
    acc_fpath = os.path.join(outdir_path, "hits_to_download.tsv")

    num_done_reads = None # variable to keep number of succeffdully processed sequences

    resume = None
    # Check if there are results from previous run.
    if os.path.exists(tsv_res_fpath) or os.path.exists(tmp_fpath):
        printl('\n' + get_work_time() + " - A result file from previous run is found in the directory:")
        printl("   '{}'".format(new_dpath))
        # Allow politely to continue from last successfully sent packet.
        resume = ask_for_resumption()
        if not resume:
            rename_file_verbosely(tsv_res_fpath, new_dpath)
            rename_file_verbosely(tmp_fpath, new_dpath)
            rename_file_verbosely(acc_fpath, new_dpath)
        # end if
    # end if
    
    # Find the name of last successfull processed sequence
    if resume == True:
        printl("Let's try to continue...")

        # Collect information from result file
        if os.path.exists(tsv_res_fpath):
            # There can be invalid information in this file
            try:
                with open(tsv_res_fpath, 'r') as res_file:
                    lines = res_file.readlines()
                    num_done_reads = len(lines) - 1 # the first line is a head
                    last_line = lines[-1]
                    last_seq_id = last_line.split(DELIM)[0]
                # end with
            except Exception as err:
                printl("\nData in classification file '{}' not found or broken. Reason:".format(tsv_res_fpath))
                printl( ' ' + str(err) )
                printl("Start from the beginning.")
                rename_file_verbosely(tsv_res_fpath, new_dpath)
                rename_file_verbosely(tmp_fpath, new_dpath)
                rename_file_verbosely(acc_fpath, new_dpath)
                return None
            else:
                printl("Last sent sequence: " + last_seq_id)
                printl("{} sequences have been already processed".format(num_done_reads))
            # end try
        # end if
        
        # Collect information from accession file
        global acc_dict
        if os.path.exists(acc_fpath):

            # There can be invalid information in this file
            try:
                with open(acc_fpath, 'r') as acc_file:
                    lines = acc_file.readlines()[9:] # omit description and head of the table
                    local_files_filtered = list( filter(lambda x: False if os.path.exists(x) else True, lines) ) # omit file paths
                    for line in local_files_filtered:
                        vals = line.split(DELIM)
                        acc = intern(vals[0].strip())
                        acc_dict[acc] = [ vals[1].strip(), vals[2].strip(), int(vals[3].strip()) ]
                    # end for
                # end with
            except Exception as err:
                printl("\nData in accession file '{}' not found or broken. Reason:".format(acc_fpath))
                printl( ' ' + str(err) )
                printl("Start from the beginning.")
                rename_file_verbosely(tsv_res_fpath, new_dpath)
                rename_file_verbosely(tmp_fpath, new_dpath)
                rename_file_verbosely(acc_fpath, new_dpath)
                return None
            else:
                printl("\nHere are Genbank records encountered during previous run:")
                for acc in acc_dict.keys():
                    s_letter = "s" if acc_dict[acc][2] > 1 else ""
                    printl(" {} hit{} - {}, '{}'".format(acc_dict[acc][2], s_letter, acc, acc_dict[acc][1]))
                # end for
                print()
            # end try
        # end if

        # If we start from the beginning, we have no sequences processed
        if num_done_reads is None:
            num_done_reads = 0
        # end if

        # Get packet size, number of the last sent packet and RID
        # There can be invalid information in tmp file of tmp file may not exist
        try:

            with open(tmp_fpath, 'r') as tmp_file:
                temp_lines = tmp_file.readlines()
            # end with

            RID_save = re_search(r"Request_ID: (.+)", temp_lines[1]).group(1).strip()
            # If aligning is performed on local machine, there is no reason for requesting results.
            # Therefore this packet will be aligned once again.
        
        except Exception as exc:

            # There is no need to disturb a user, merely resume.
            return {
                "RID": None,
                "acc_fpath": acc_fpath,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        else:
            # Return data from previous run
            return {
                "RID": RID_save,
                "acc_fpath": acc_fpath,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        # end try
    
    return None
# end def look_around

# According to
# https://github.com/nanoporetech/ont_h5_validator/blob/master/h5_validator/schemas/multi_read_fast5.yaml
ont_read_signature = r"([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})"

def fmt_read_id(read_id):

    srch_ont_read = re_search(ont_read_signature, read_id)
    if srch_ont_read is None:
        return read_id.partition(' ')[0]
    else:
        return '>' + srch_ont_read.group(1)
# end def fmt_read_id


def pass_processed_seqs(fasta_file, num_done_reads, fmt_func):
    """
    Function passes sequences that have been already processed.

    :param fasta_file: FASTA file instalce;
    :type fasta_file: str;
    :param num_done_reads: amount of sequences that have been already processed;
    :type num_done_reads: int;
    :param fmt_func: function from 'FORMATTING_FUNCS' tuple;
    """

    if num_done_reads == 0:
        return None
    else:
        i = 0
        while i <= num_done_reads:

            line = fmt_func(fasta_file.readline())
            if line == "":
                return ""
            # end if
            if line.startswith('>'):
                line = fmt_read_id(line)
                next_id_line = line
                i += 1
            # end if
        # end while
        return next_id_line
    # end if
# end def pass_processed_seqs


def fasta_packets(fasta, packet_size, reads_at_all, num_done_reads):
    """
    Function (actually, generator) retrieves 'packet_size' records from FASTA file.
    This function will pass 'num_done_reads' sequences (i.e. they will not be processed)
        by calling 'pass_processed_files'.

    :param fasta: path to FASTA file;
    :type fasta: str;
    :param packet_size: number of sequences to align in one 'blastn' launching;
    :type packet_size: int;
    :param reads_at_all: number of sequences in current file;
    :type reads_at_all: int;
    :param num_done_reads: number of sequnces in current file that have been already processed;
    :type num_doce_reads: int;
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fasta) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fasta) ]

    with how_to_open(fasta) as fasta_file:
        # Next line etrieving will be performed as simple line-from-file reading.
        get_next_line = lambda: fmt_func(fasta_file.readline())

        # Variable that contains id of next sequence in current FASTA file.
        # If no or all sequences in current FASTA file have been already processed, this variable is None.
        # There is no way to count sequences in multi-FASTA file, accept of counting sequence IDs.
        # Therefore 'next_id_line' should be saved in memory after moment when packet is formed.
        try:
            next_id_line = pass_processed_seqs(fasta_file, num_done_reads, fmt_func)
        except UnboundLocalError:
            # This exception occurs when 'fasta_file' variable is not defined, i.e. when
            #   'fasta' is actual FASTA data, not path to file.
            # In this case we need all FASTA data.
            next_id_line = None
        # end try

        if next_id_line == "":
            yield {"fasta": "", "names": list()}
        # end if

        packet = ""

        line = get_next_line()
        if line.startswith('>'):
            line = fmt_read_id(line) # prune sequence ID
        # end if

        # If some sequences have been passed, this if-statement will be executed.
        # New packet should start with sequence ID line.
        if not next_id_line is None:
            packet += next_id_line+'\n'
        # end if
        packet += line+'\n' # add recently read line

        packs_at_all = reads_at_all // packet_size # Calculate total number of packets sent from current FASTA file
        if reads_at_all % packet_size > 0: # And this is ceiling (in order not to import 'math')
            packs_at_all += 1
        # end if
        packs_processed = int( num_done_reads / packet_size ) # number of successfully processed sequences
        packs_left = packs_at_all - packs_processed # number of packets left to send

        # Iterate over packets left to process
        for _ in range(packs_left):

            i = 0 # variable for counting sequenes within packet

            packet_size = min(packet_size, probing_batch_size - seqs_processed)
            
            while i < packet_size:

                line = get_next_line()
                if line.startswith('>'):
                    line = fmt_read_id(line)
                    i += 1
                # end if
                
                if line == "": # if end of file (data) is reached
                    break
                # end if
                packet += line + '\n' # add line to packet
            # end while

            if line != "":
                next_id_line = packet.splitlines()[-1] # save sequence ID next packet will start with
                packet = '\n'.join(packet.splitlines()[:-1]) # exclude 'next_id_line' from packet
            else:
                next_id_line = None
            # end if

            # Get list of sequence IDs:
            names = list( filter(lambda l: True if l.startswith('>') else False, packet.splitlines()) )
            names = list( map(lambda l: l.replace('>', ''), names) )

            # Just in case
            if packet == "":
                printl("Recent packet is empty")
                return
            # end if

            if max_seq_len is None:
                yield {"fasta": packet.strip(), "names": names}
            else:
                yield {"fasta": prune_seqs(packet.strip(), 'l', max_seq_len),
                    "names": names}
            # end if

            # Reset packet
            if not next_id_line is None:
                packet = next_id_line+'\n'
            else:
                return
            # end if
        # end for
    # end with
# end def fasta_packets


def remove_tmp_files(*paths):
    """
    Function removes files passed to it.
    Actually, passed arguments are paths ('str') to files meant to be removed.
    """

    for path in paths:
        if os.path.exists(path):
            os.unlink(path)
        # end if
    # end for
# end def remove_tmp_files


def configure_request(packet, blast_algorithm, organisms):
    """
    Function configures the request to BLAST server.

    :param packet: FASTA_data_containing_query_sequences;
    :type packet: str;
    :param blast_algorithm: BLASTn algorithm to use;
    :type blast_algorithm: str;
    :param organisms: list of strings performing 'nt' database slices;
    :type organisms: list<str>;

    Returns a dict of the following structure:
    {
        "payload": the_payload_of_the_request (dict),
        "headers": headers of thee request (dict)
    }
    """

    payload = dict()
    payload["CMD"] = "PUT" # method
    payload["PROGRAM"] = "blastn" # program
    payload["MEGABLAST"] = "on" if "megablast" in blast_algorithm.lower() else "" # if megablast
    payload["BLAST_PROGRAMS"] = blast_algorithm # blastn algoeithm
    payload["DATABASE"] = "nt" # db
    payload["QUERY"] = packet # FASTA data
    payload["HITLIST_SIZE"] = 1 # we need only the best hit
    if user_email != "":
        payload["email"] = user_email # user's email
        payload["tool"] = "barapost:_prober"
    # end if
    

    # 'nt' database slices:
    for i, org in enumerate(organisms):
        payload["EQ_MENU{}".format(i if i > 0 else "")] = org
    # end for
    
    payload["NUM_ORG"] = str( len(organisms) )

    payload = urllib.parse.urlencode(payload)

    headers = { "Content-Type" : "application/x-www-form-urlencoded" }

    return {"payload":payload, "headers": headers}
# end def configure_request


def send_request(request, pack_to_send, packs_at_all, filename, tmp_fpath):
    """
    Function sends a request to "blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        and then waits for satisfaction of the request and retrieves response text.

    :param request: request_data (it is a dict that 'configure_request()' function returns);
    :param request: dict<dict>;
    :param pack_to_send: current number (like id) of packet meant to be sent now.
    :type pack_to_send: int;
    :param packs_at all: total number of packets corresponding to current FASTA file.
        This information is printed to console;
    :type packs_at_all: int;
    :param filename: basename of current FASTA file;
    :type filename: str;

    Returns XML text of type 'str' with BLAST response.
    """
    payload = request["payload"]
    headers = request["headers"]

    server = "blast.ncbi.nlm.nih.gov"
    url = "/blast/Blast.cgi"
    error = True

    while error:
        try:
            conn = http.client.HTTPSConnection(server) # create a connection
            conn.request("POST", url, payload, headers) # send the request
            response = conn.getresponse() # get the response
            response_text = str(response.read(), "utf-8") # get response text
        except OSError as oserr:
            printl(get_work_time() + "\n - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            printl( repr(err) )
            printl("barapost will try to connect again in 30 seconds...\n")
            sleep(30)
        
        # if no exception occured
        else:
            error = False
        # end try
    # end while

    try:
        rid = re_search(r"RID = (.+)", response_text).group(1) # get Request ID
        rtoe = int(re_search(r"RTOE = ([0-9]+)", response_text).group(1)) # get time to wait provided by the NCBI server
    except AttributeError:
        printl(err_fmt("seems, ncbi has denied your request."))
        printl("Response is in file 'request_denial_response.html'")
        with open("request_denial_response.html", 'w') as den_file:
            den_file.write(response_text)
        # end with
        platf_depend_exit(1)
    finally:
        conn.close()
    # end try

    # Save temporary data
    with open(tmp_fpath, 'w') as tmpfile:
        tmpfile.write("sent_packet_num: {}\n".format(pack_to_send))
        tmpfile.write("Request_ID: {}".format(rid))
    # end with

    # /=== Wait for alignment results ===\

    return( wait_for_align(rid, rtoe, pack_to_send, packs_at_all, filename) )
# end def send_request


def lingering_https_get_request(server, url):
    """
    Function performs a "lingering" HTTPS request.
    It means that the function tries to get the response
        again and again if the request fails.

    :param server: server address;
    :type server: str;
    :param url: the rest of url;
    :type url: str;

    Returns obtained response decoded to UTF-8 ('str').
    """

    error = True
    while error:
        try:
            conn = http.client.HTTPSConnection(server) # create connection
            conn.request("GET", url) # ask for if there areresults
            response = conn.getresponse() # get the resonse
            resp_content = str(response.read(), "utf-8") # get response text
            conn.close()
        except OSError as err:
            printl("\n{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
            printl("  " + str(err))
            error = True
            sleep(30)
        except http.client.RemoteDisconnected as err:
            printl("\n{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
            printl("  " + str(err))
            error = True
            sleep(30)
        except socket.gaierror as err:
            printl("\n{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
            printl("  " + str(err))
            error = True
            sleep(30)
        except http.client.CannotSendRequest as err:
            printl("\n{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
            printl("  " + str(err))
            error = True
            sleep(30)
        else:
            error = False # if no exception ocured
        # end try
    # end while
    return resp_content

# end def lingering_https_get_request


def wait_for_align(rid, rtoe, pack_to_send, packs_at_all, filename):
    """
    Function waits untill BLAST server accomplishes the request.
    
    :param rid: Request ID to wait for;
    :type rid: str;
    :param rtoe: time in seconds estimated by BLAST server needed to accomplish the request;
    :type rtoe: int;
    :param pack_to_send: current packet (id) number to send;
    :type pack_to_send: int;
    :param packs_at_all: total number of packets corresponding to current FASTA file.
        This information is printed to console;
    :type packs_at_all: int;
    :param filename: basename of current FASTA file;
    :type filename: str;

    Returns XML response ('str').
    """

    printl("\n{} - Requesting for alignment results. Request ID: {},\n '{}' ({}/{})".format(get_work_time(),
    rid, filename, pack_to_send, packs_at_all))
    # RTOE can be zero at the very beginning of resumption
    if rtoe > 0:
        printl("{} - BLAST server estimates that alignment will be accomplished in {} seconds ".format(get_work_time(), rtoe))
        printl("{} - Waiting for {}+3 (+3 extra) seconds...".format(get_work_time(), rtoe))
        # Server migth be wrong -- we will give it 3 extra seconds
        sleep(rtoe + 3)
        printl("{} - {} seconds have passed. Checking if alignment is accomplished...".format(get_work_time(), rtoe+3))
    # end if

    server = "blast.ncbi.nlm.nih.gov"
    wait_url = "/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid
    there_are_hits = False

    while True:
        resp_content = lingering_https_get_request(server, wait_url)

        # if server asks to wait
        if "Status=WAITING" in resp_content:
            printn("\r{} - The request is still processing. Waiting      \033[6D".format(get_work_time()))
            # indicate each 20 seconds with a dot
            for i in range(1, 7):
                sleep(10)
                printn("\r{} - The request is still processing. Waiting{}".format(get_work_time(), '.'*i))
            # end for
            continue
        # end if
        if "Status=FAILED" in resp_content:
            # if job failed
            printl('\n' + get_work_time() + " - Job failed\a\n")
            return """{} - Job for query {} ({}/{}) with Request ID {} failed.
    Contact NCBI or try to start it again.\n""".format(get_work_time(), filename, pack_to_send, packs_at_all, rid)
        # end if
        # if job expired
        if "Status=UNKNOWN" in resp_content:
            printl('\n' + get_work_time() + " - Job expired\a\n")
            return "expired"
        # end if
        # if results are ready
        if "Status=READY" in resp_content:
            there_are_hits = True
            printl("\n{} - Result for query '{}' ({}/{}) is ready!".format(get_work_time(), filename, pack_to_send, packs_at_all))
            # if there are hits
            if "ThereAreHits=yes" in resp_content:
                for i in range(45, 0, -5):
                    printl('-' * i)
                # end for
                print() # just print a blank line
            # if there are no hits
            else:
                printl(get_work_time() + " - There are no hits. It happens.\n")
            # end if
            break
        # end if
        # Execution should not reach here
        printl('\n' + get_work_time() + " - Fatal error. Please contact the developer.\a\n")
        platf_depend_exit(1)
    # end while

    # Retrieve XML result
    retrieve_xml_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&ALIGNMENTS=1&RID=" + rid
    respond_text= lingering_https_get_request(server, retrieve_xml_url)

    if "[blastsrv4.REAL]" in respond_text:
        printl("BLAST server error:\n  {}".format(re_search(r"(\[blastsrv4\.REAL\].*\))", respond_text).group(1)))
        printl("All sequences in recent packet will be halved and this packet will be resent to NCBI BLAST server.")
        return None
    # end if

    # Retrieve human-readable text and put it into result directory
    if there_are_hits:
        save_txt_align_result(server, filename, pack_to_send, rid)
    # end if

    return respond_text
# end def wait_for_align


def prune_seqs(packet, mode, value):
    """
    Function prunes all sequences in packet, leaving 5'-half of the sequence.

    :param packet: FASTA data;
    :type packet: str;
    :param mode: available values: 'l' (for length) and 'f' (for fraction).
        Meaning: if mode is 'l', all sequences are pruned from the 5'-end to position 'value' (value > 1),
                 if mode is 'f', all sequences are pruned from the 5'-end to position L*value (0<value<1),
    :type mode: str;
    :param value: see description if 'mode' parameter above;
    :type value:str;
    """

    lines = packet.splitlines()
    packet = ""

    if mode == 'l':
        if value < 1:
            raise ValueError("Invalid value of parameter 'value': it must be > 1")
        # end if

        def prune(seq, value):
            try:
                return seq[: value]
            except IndexError:
                return seq
            # end try
        # end def prune
    elif mode == 'f':
        if value < 0.0 + 1e-9 or value > 1.0 - 1e-9:
            raise ValueError("Invalid value of parameter 'value': it must be in interval (0, 1)")
        # end if

        def prune(seq, value):
            return seq[: int(len(seq) * value)]
        # end def prune
    else:
        raise ValueError("Invalid 'mode' argument passed to function prune_seqs: '{}'".format(mode))
    # end if

    id_line_idxs = list(map(lines.index, filter(lambda x: x.startswith('>'), lines)))
    id_line_idxs.append(len(lines)) # fictive last id line

    for id_curr, id_next in zip(id_line_idxs[:-1], id_line_idxs[1:]):
        is_seq = lambda seq_i: seq_i > id_curr and seq_i < id_next
        seq_lines = map(lambda j: lines[j], filter(is_seq, range(len(lines))))
        seq = "".join(seq_lines)
        seq = prune(seq, value)
        packet += lines[id_curr] + '\n' + seq + '\n'
    # end for

    return packet
# end def prune_seqs


def save_txt_align_result(server, filename, pack_to_send, rid):

    global outdir_path

    retrieve_text_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&DESCRIPTIONS=1&ALIGNMENTS=1&RID=" + rid
    respond_text = lingering_https_get_request(server, retrieve_text_url)

    txt_hpath = os.path.join(outdir_path, "prober_blast_response_{}.txt".format(pack_to_send))
    # Write text result for a human to read
    with open(txt_hpath, 'w') as txt_file:
        txt_file.write(respond_text + '\n')
    # end with

# end def save_txt_align_result


def get_lineage(gi, hit_def, hit_acc):
    """
    Function retrieves lineage of a hit from NCBI.
    It downloads INSDSeq XML file, since it is the smallest one among those containing lineage.
    Moreover, it saves this lineage in 'taxonomy' DBM file:
        {<accession>: <lineage_str>}

    :param gi: GI number of a hit;
    :type gi: str;
    :param hit_def: definition line of a hit;
    :type hit_def: str;
    :param hit_acc: hit accession;
    :type hit_acc: str;
    """

    # Get all accessions in taxonomy file:
    tax_acc_exist = shelve.open(taxonomy_path, 'c').keys()

    # If we've got a new accession -- download lineage
    if not hit_acc in tax_acc_exist:

        retrieve_url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbc_xml&id={}&".format(gi)

        error = True
        while error:
            try:
                urllib.request.urlretrieve(retrieve_url, indsxml_path)
            except OSError:
                print("\nError while requesting for lineage.\n Let's try again in 30 seconds.")
                if os.path.exists(indsxml_path):
                    os.unlink(indsxml_path)
                # end if
                sleep(30)
            else:
                error = False
            # end try
        # end while

        # Get downloaded text:
        text = open(indsxml_path, 'r').read()

        try:
            # Find genus name and species name:
            org_name_regex = r"<INSDSeq_organism>([A-Z][a-z]+ [a-z]+(\. [a-zA-Z0-9]+)?)"
            org_name = re_search(org_name_regex, text).group(1)

            # Get species name:
            spec_name = re_search(r" ([a-z]+(\. [a-zA-Z0-9]+)?)", org_name).group(1)

            # Get full lineage
            lin_regex = r"<INSDSeq_taxonomy>([A-Za-z; ]+)</INSDSeq_taxonomy>"
            lineage = re_search(lin_regex, text).group(1).strip('.')

            # Check if species name is in lineage:
            gen_spec_regex = r"([A-Z][a-z]+ [a-z]+).?$"
            gen_spec_in_lin = re_search(gen_spec_regex, lineage)

            # All taxonomy names in lineage will be divided by semicolon:
            if not gen_spec_in_lin is None:
                repl_str = gen_spec_in_lin.group(1)
                lineage = lineage.replace(repl_str, spec_name)
            else:
                lineage += ";" + spec_name
            # emd if

            # Remove all spaces:
            lineage = lineage.replace(' ', '')

        except AttributeError:
            # If there is no lineage -- use hit definition instead of it

            # Format hit definition (get rid of stuff after comma)
            hit_def = hit_def[: hit_def.find(',')] if ',' in hit_def else hit_def
            hit_def = hit_def.replace(" complete genome", "") # sometimes there are no comma before it
            hit_def = hit_def.replace(' ', '_')

            shelve.open(taxonomy_path, 'c')[hit_acc] = hit_def
            lineage = hit_def
        # end try

        # Write lineage to taxonomy file
        with shelve.open(taxonomy_path, 'c') as tax_file:
            tax_file[str(hit_acc)] = lineage
    else:
        # If hit is not new -- simply retrieve it from taxonomy file
        with shelve.open(taxonomy_path, 'c') as tax_file:
            lineage = tax_file[str(hit_acc)]
    # end if

    # Remove tmp xml file
    if os.path.exists(indsxml_path):
        os.unlink(indsxml_path)
    # end if

    return lineage
# end def get_lineage


def parse_align_results_xml(xml_text, seq_names):
    """
    Function parses BLAST xml response and returns tsv lines containing gathered information:
        1. Query name.
        2. Hit name formatted by 'format_taxonomy_name()' function.
        3. Hit accession.
        4. Length of query sequence.
        5. Length of alignment.
        6. Percent of identity.
        7. Percent of gaps.
        8. E-value.
        9. Average Phred33 quality of a read (if source file is FASTQ).
        10. Read accuracy (%) (if source file is FASTQ).

    Erroneous tsv lines that function may produce:
        1. "<query_name>\\tQuery has been lost: ERROR, Bad Gateway"
            if data packet has been lost.
            # end if
        2. "<query_name>\\tQuery has been lost: BLAST ERROR"
            if BLAST error occured.
            # end if
        3. "<query_name>\\tNo significant similarity found"
            if no significant similarity has been found
            # end if
        Type of return object: list<str>.
    """

    result_tsv_lines = list()

    # /=== Validation ===/

    if "Bad Gateway" in xml_text:
        printl('\n' + '=' * 45)
        printl(get_work_time() + " - ERROR! Bad Gateway! Data from last packet has lost.")
        printl("It would be better if you restart the script.")
        printl("Here are names of lost queries:")
        for i, name in enumerate(seq_names):
            printl("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + DELIM + "Query has been lost: ERROR, Bad Gateway")
        # end for
        return result_tsv_lines
    # end if

    if "Status=FAILED" in xml_text:
        printl('\n' + get_work_time() + "BLAST ERROR!: request failed")
        printl("Here are names of lost queries:")
        for i, name in enumerate(seq_names):
            printl("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + DELIM +"Query has been lost: BLAST ERROR")
        # end for
        return result_tsv_lines
    # end if

    if "to start it again" in xml_text:
        printl('\n' + get_work_time() + "BLAST ERROR!")

        printl("Here are names of lost queries:")
        for i, name in enumerate(seq_names):
            printl("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + DELIM +"Query has been lost: BLAST ERROR")
        # end for
        return result_tsv_lines
    # end if

    # /=== Parse BLAST XML response ===/

    root = ElementTree.fromstring(xml_text) # get tree instance

    global new_acc_dict
    global qual_dict

    # Iterate through "Iteration" and "Iteration_hits" nodes
    for iter_elem, iter_hit in zip(root.iter("Iteration"), root.iter("Iteration_hits")):
        # "Iteration" node contains query name information
        query_name = intern(iter_elem.find("Iteration_query-def").text)
        query_len = iter_elem.find("Iteration_query-len").text

        if not qual_dict is None:
            ph33_qual = qual_dict['>' + query_name]
            miscall_prop = round(10**(ph33_qual/-10), 3)
            accuracy = round( 100*(1 - miscall_prop), 2 ) # expected percent of correctly called bases
            qual_info_to_print = "    Average quality of this read is {}, i.e. accuracy is {}%;\n".format(ph33_qual,
                accuracy)
        else:
            # If FASTA file is processing, print dashed in quality columns
            ph33_qual = "-"
            accuracy = "-" # expected percent of correctly called bases
            qual_info_to_print = ""
        # end if

        # Check if there are any hits
        chck_h = iter_hit.find("Hit")

        if chck_h is None:
            # If there is no hit for current sequence
            printl("\n '{}' -- No significant similarity found;\n    Query length - {};".format(query_name, query_len))
            result_tsv_lines.append(DELIM.join( (query_name, "No significant similarity found", "-", query_len,
                "-", "-", "-", "-", str(ph33_qual), str(accuracy)) ))
        else:
            # If there are any hits, node "Iteration_hits" contains at least one "Hit" child
            # Get first-best bitscore and iterato over hits that have the save (i.e. the highest bitscore):
            top_bitscore = next(chck_h.find("Hit_hsps").iter("Hsp")).find("Hsp_bit-score").text

            lineages = list()
            hit_accs = list()

            for hit in iter_hit:

                # Find the first HSP (we need only the first one)
                hsp = next(hit.find("Hit_hsps").iter("Hsp"))

                if hsp.find("Hsp_bit-score").text != top_bitscore:
                    break
                # end if

                # Get full hit name (e.g. "Erwinia amylovora strain S59/5, complete genome")
                hit_def = hit.find("Hit_def").text

                curr_acc = intern(hit.find("Hit_accession").text)
                hit_accs.append( curr_acc ) # get hit accession
                gi_patt = r"gi\|([0-9]+)" # pattern for GI number finding
                hit_gi = re_search(gi_patt, hit.find("Hit_id").text).group(1)

                # Get lineage
                try:
                    lineages.append(get_lineage(hit_gi, hit_def, hit_accs[len(hit_accs)-1]).split(';'))
                except OSError as oserr:
                    print(err_fmt(str(oserr)))
                    platf_depend_exit(1)
                # end try

                # Update accession dictionary
                try:
                    new_acc_dict[curr_acc][2] += 1
                except KeyError:
                    new_acc_dict[curr_acc] = [hit_gi, hit_def, 1]
                # end try

                align_len = hsp.find("Hsp_align-len").text.strip()
                pident = hsp.find("Hsp_identity").text # get number of matched nucleotides
                gaps = hsp.find("Hsp_gaps").text # get number of gaps

                evalue = hsp.find("Hsp_evalue").text # get e-value
                pident_ratio = round( float(pident) / int(align_len) * 100, 2)
                gaps_ratio = round( float(gaps) / int(align_len) * 100, 2)
            # end for

            # Finc LCA:
            lineage = list()

            for i in range(len(lineages[0])):
                if len(set(map(lambda t: t[i], lineages))) == 1:
                    lineage.append(lineages[0][i]) 
                else:
                    break
                # end if
            # end for

            if len(lineage) != 0:
                # Divide taxonomic names with semicolon
                lineage = ';'.join(lineage)
                printl("""\n '{}' -- '{}';
    Query length - {} nt;
    Identity - {}/{} ({}%); Gaps - {}/{} ({}%);""".format(query_name, lineage,
                    query_len, pident, align_len, pident_ratio, gaps, align_len, gaps_ratio))
                # Append new tsv line containing recently collected information
                result_tsv_lines.append( DELIM.join( (query_name, lineage, '&&'.join(hit_accs), query_len,
                    align_len, pident, gaps, evalue, str(ph33_qual), str(accuracy)) ))
            else:
                # If hits are from different domains (i.e. Bacteria, Archaea, Eukaryota) -- no let it be unknown
                printl("\n '{}' -- No significant similarity found;\n    Query length - {};".format(query_name, query_len))
                result_tsv_lines.append(DELIM.join( (query_name, "No significant similarity found", "-", query_len,
                    "-", "-", "-", "-", str(ph33_qual), str(accuracy)) ))
            # end if
        # end if
        println(qual_info_to_print)
    # end for

    return result_tsv_lines
# end def parse_align_results_xml


def write_result(res_tsv_lines, tsv_res_path, acc_file_path, fasta_hname, outdir_path):
    """
    Function writes result of blasting to result tsv file.

    :param res_tsv_lines: tsv lines returned by 'parse_align_results_xml()' funciton;
    :type res_tsv_lines: list<str>;
    :param tsv_res_path: path to reslut tsv file;
    :type tsv_res_path: str;
    :param new_dpath: path to current result directory;
    :type new_dpath: str;
    """

    # If there is no result tsv file -- create it and write a head of the table.
    if not os.path.exists(tsv_res_path):
        with open(tsv_res_path, 'w') as tsv_res_file:
            tsv_res_file.write(DELIM.join( ["QUERY_ID", "HIT_NAME", "HIT_ACCESSION", "QUERY_LENGTH",
                "ALIGNMENET_LENGTH", "IDENTITY", "GAPS", "E-VALUE", "AVG_PHRED33", "ACCURACY(%)"] ) + '\n')
        # end with
    # end if
    
    # Write reslut tsv lines to this file
    with open(tsv_res_path, 'a') as tsv_res_file:
        for line in res_tsv_lines:
            tsv_res_file.write(line + '\n')
        # end for
    # end with

    # === Write accession information ===

    global acc_dict
    global new_acc_dict
    global blast_algorithm
    acc_file_path = os.path.join(outdir_path, "hits_to_download.tsv")

    with open(acc_file_path, 'w') as acc_file:
        acc_file.write("# Here are accessions, GI numbers and descriptions of Genbank records that can be used for sorting by 'barapost.py'\n")
        acc_file.write("# Values in this file are delimited by tabs.\n")
        acc_file.write("# You are welcome to edit this file by adding,\n")
        acc_file.write("#   removing or muting lines (with adding '#' symbol in it's beginning, just like this description).\n")
        acc_file.write("# Lines muted with '#' won't be noticed by 'barapost.py'.\n")
        acc_file.write("# You can specify your own FASTA files that you want to use as database for 'barapost.py'.\n")
        acc_file.write("# To do it, just write your FASTA file's path to this TSV file in new line.\n\n")
        acc_file.write(DELIM.join( ["ACCESSION", "GI_NUMBER", "RECORD_NAME", "OCCURRENCE_NUMBER"] ) + '\n')
    # end with

    for acc, other_info in new_acc_dict.items():
        try:
            acc_dict[acc][2] += other_info[2]
        
        except KeyError:
            acc_dict[acc] = other_info
        # end try
    # end for
    
    # Write accessions and record names
    with open(acc_file_path, 'a') as acc_file:
        for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][2]):
            acc_file.write(DELIM.join( (acc, other_info[0], other_info[1], str(other_info[2]))) + '\n')
        # end for
    # end with
    
    new_acc_dict.clear()
# end def write_result


def create_result_directory(fq_fa_path, outdir_path):
    """
    Function creates a result directory named according 
        to how source FASTQor FASTA file is named.

    :param fq_fa_path: path to source FASTQ or FASTA file;
    :type fq_fa_path: str;

    Returns 'str' path to the recently created result directory.
    """

    # dpath means "directory path"
    new_dpath = os.path.join(outdir_path, os.path.basename(fq_fa_path)) # get rid of absolute path
    new_dpath =re_search(r"(.*)\.(m)?f(ast)?(a|q)", new_dpath).group(1) # get rid of extention
    if not os.path.exists(new_dpath):
        try:
            os.makedirs(new_dpath)
        
        except OSError as oserr:
            printl(err_fmt("error while creating result directory"))
            printl( str(oserr) )
        # end try
    # end if

    return new_dpath
# end def create_result_directory


# =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=
#                       |===== Proceed =====|


# /=== Comments to the kernel loop ===/

# 1. 'curr_fasta' is a dict of the following structure:
#    {
#        "fpath": path_to_FASTA_file (str),
#        "nreads": number_of_reads_in_this_FASTA_file (int)
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. 'previous_data' is a dict of the following structure:
#    {
#        "sv_npck": saved_number_of_sent_packet (int),
#        "RID": saved_RID (str),
#        "tsv_respath": path_to_tsv_file_from_previous_run (str),
#        "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file (int),
#        "tmp_fpath": path_to_pemporary_file (str)
#        "acc_fparh": path_to_accessoion_file (str)
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. 'packet' is a dict of the following structure:
#    {
#        "fasta": FASTA_data_containing_query_sequences (str),
#        "names": list_of_sequence_ids_from_FASTA_file (list<str>)
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. 'response' is a dict of the following structure:
#    {
#        "RID": Request ID (str),
#        "RTOE", time_to_wait_provided_by_the_NCBI_server (int)
#    }

#                   |===== Kernel loop =====|

organisms = verify_taxids(taxid_list)

printl(" - Output directory: '{}';".format(outdir_path))
if user_email != "":
    printl(" - Your email: {}".format(user_email))
# end if
printl(" - Probing batch size: {} sequences;".format(probing_batch_size))
printl(" - Packet size: {} sequences;".format(packet_size))
printl(" - BLAST algorithm: {};".format(blast_algorithm))
if not max_seq_len is None:
    printl(" - Maximum length of a sequence sent: {} bp;".format(max_seq_len))
# end if
printl(" - Database: nt;")
if len(organisms) > 0:
    for db_slice in organisms:
        printl("   {};".format(db_slice))
    # end for
# end if

s_letter = '' if len(fq_fa_list) == 1 else 's'
printl("\n {} file{} will be processed.".format( len(fq_fa_list), s_letter))
logfile.write("Here they are:\n")
for i, path in enumerate(fq_fa_list):
    logfile.write("    {}. '{}'\n".format(i+1, path))
# end for

printl('-'*30)

# Variable for counting accessions of records meant to be downloaded from Genbank.
# Is used only for printing the list of accessions to console.
acc_counter = 0
# Dictionary of accessions and record names.
# Accessions are keys, record names are values.
# This dictionary is filled while processing and at the beginning of resumption.
acc_dict = dict()
# Dictionary of accessions and record names encountered while sending of the current packet.
# Accessions are keys, record names are values.
new_acc_dict = dict()

# Counter of sequences processed during current run
seqs_processed = 0
# Counetr of sequences processed concerning putative previous run(s)
glob_seqs_processed = 0

# Variable that contains id of next sequence in current FASTA file.
# If no or all sequences in current FASTA file have been already processed, this variable is None
# Function 'get_packet' changes this variable
next_id_line = None

# Varible for stopping execution when probing batch is processed completely.
stop = False

# Iterate through found source FASTQ and FASTA files
for i, fq_fa_path in enumerate(fq_fa_list):

    # Configure quality dictionary
    qual_dict = configure_qual_dict(fq_fa_path) if is_fastq(fq_fa_path) else None

    # Create the result directory with the name of FASTQ of FASTA file being processed:
    new_dpath = create_result_directory(fq_fa_path, outdir_path)

    # Convert FASTQ file to FASTA (if it is FASTQ) and get it's path and number of sequences in it:
    curr_fasta = fastq2fasta(fq_fa_path, i, new_dpath)
    printl("\n |=== file: '{}' ({} sequences) ===|".format(os.path.basename(curr_fasta["fpath"]),
        curr_fasta["nreads"]))

    # "hname" means human readable name (i.e. without file path and extention)
    fasta_hname = os.path.basename(curr_fasta["fpath"]) # get rid of absolure path
    fasta_hname = re_search(r"(.*)\.(m)?f(ast)?a", fasta_hname).group(1) # get rid of file extention

    # Look around and ckeck if there are results of previous runs of this script
    # If 'look_around' is None -- there is no data from previous run
    previous_data = look_around(outdir_path, new_dpath, curr_fasta["fpath"],
        blast_algorithm)

    if previous_data is None: # If there is no data from previous run
        num_done_reads = 0 # number of successfully processed sequences
        tsv_res_path = "{}.tsv".format(os.path.join(new_dpath,
            "classification")) # form result tsv file path
        tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form temporary file path
        acc_fpath = os.path.join(outdir_path, "hits_to_download.tsv") # form path to accession file
        saved_RID = None
    else: # if there is data from previous run
        num_done_reads = previous_data["n_done_reads"] # get number of successfully processed sequences
        tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        tmp_fpath = previous_data["tmp_fpath"] # temporary file sholud be the same as during previous run
        acc_fpath = previous_data["acc_fpath"] # accession file sholud be the same as during previous run
        saved_RID = previous_data["RID"] # having this RID we can try to get response for last request
        contin_rtoe = 0 # we will not sleep at the very beginning of resumption
    # end if

    glob_seqs_processed += num_done_reads

    # Calculate total number of packets meant to be sent from current FASTA file
    tmp_num = min(probing_batch_size, curr_fasta["nreads"])
    packs_at_all = tmp_num // packet_size

    if tmp_num % packet_size > 0: # And this is ceiling (in order not to import 'math')
        packs_at_all += 1
    # end if

    pack_to_send = 1 # number of packet meant to be sent now

    if is_gzipped(curr_fasta["fpath"]):
        fmt_func = lambda l: l.decode("utf-8")
        how_to_open = open_as_gzip
    
    else:
        fmt_func = lambda l: l
        how_to_open = open
    # end if

    # Iterate over packets left to send
    for packet in fasta_packets(curr_fasta["fpath"], packet_size, curr_fasta["nreads"], num_done_reads):

        # Just in case:
        if packet["fasta"] == "":
            printl("All sequences in recent file have been already processed.")
            break
        # end if

        send = True

        # If current packet has been already send, we can try to get it and not to send it again
        if not saved_RID is None:

            align_xml_text = wait_for_align(saved_RID, contin_rtoe,
                pack_to_send, packs_at_all, fasta_hname+".fasta") # get BLAST XML response
            saved_RID = None

            # If request is not expired get he result and not send it again
            if align_xml_text != "expired":
                send = False

                result_tsv_lines = parse_align_results_xml(align_xml_text,
                    packet["names"]) # get result tsv lines

                seqs_processed += len( packet["names"] )

                # Write the result to tsv
                write_result(result_tsv_lines, tsv_res_path, acc_fpath, fasta_hname, outdir_path)
            # end if
        # end if

        if send:

            printl("\nGo to BLAST (" + blast_algorithm + ")!")
            printl("Request number {} out of {}. Sending {} sequences.".format(pack_to_send,
                packs_at_all, len(packet["names"])))

            align_xml_text = None
            while align_xml_text is None: # untill successfull attempt

                request = configure_request(packet["fasta"], blast_algorithm, organisms) # get the request

                # Send the request get BLAST XML response
                # 'align_xml_text' will be None if NCBI BLAST server rejects the request due to too large amount of data in it.

                align_xml_text = send_request(request,
                    pack_to_send, packs_at_all, fasta_hname+".fasta", tmp_fpath)

                if align_xml_text == "expired":
                    printl("Job expired. Trying to send it once again.\n")
                    continue
                # end if

                # If NCBI BLAST server rejects the request due to too large amount of data in it --
                #    shorten all sequences in packet twofold and resend it.
                if align_xml_text is None:
                    packet["fasta"] = prune_seqs(packet["fasta"], 'f', 0.5)
                # end if
            # end while

            # Get result tsv lines
            result_tsv_lines = parse_align_results_xml(align_xml_text,
                packet["names"])

            seqs_processed += len( packet["names"] )

            # Write the result to tsv
            write_result(result_tsv_lines, tsv_res_path, acc_fpath, fasta_hname, outdir_path)
        # end if
        pack_to_send += 1
        remove_tmp_files(tmp_fpath)

        if seqs_processed >= probing_batch_size:
            stop = True
            break
        # end if
    # end for
    remove_tmp_files(tmp_fpath)
    if stop:
        break
    # end if
# end for

if os.path.exists(indsxml_path):
    os.unlink(indsxml_path)
# end if

def get_undr_sep_number(number):
    undr_sep_num = str(number)
    for i in range(len(undr_sep_num)-4, -1, -4):
        undr_sep_num = undr_sep_num[: i+1] + '_' + undr_sep_num[i+1: ]
    # end for
    return undr_sep_num
# end def get_undr_sep_number

glob_seqs_processed += seqs_processed
str_about_prev_runs = ", including previous run(s)" if glob_seqs_processed > seqs_processed else ""

printl('-'*20+'\n')
printl(" {} sequences have been processed{}\n".format(get_undr_sep_number(glob_seqs_processed),
    str_about_prev_runs))

printl("Here are Genbank records that can be used for further sorting by 'barapost.py'.")
printl("They are sorted by their occurence in probing batch:")

# Print accessions and record names sorted by occurence
# "-x[1][2]:": minus because we need descending order, [1] -- get tuple of "other information",
#   [2] -- get 3-rd element (occurence)
for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][2]):
    s_letter = "s" if other_info[2] > 1 else ""
    printl(" {} hit{} - {}, '{}'".format(get_undr_sep_number(other_info[2]),
        s_letter, acc, other_info[1]))
# end for

# Print number of unkmown sequences, if there are any:
unkn_num = glob_seqs_processed - sum( map(lambda x: x[2], acc_dict.values()) )
if unkn_num > 0:
    s_letter = "s" if unkn_num > 1 else ""
    printl(" {} hit{} - No significant similarity found".format(unkn_num, s_letter))
# end if

printl("""\nThey are saved in following file:
    '{}'""".format(acc_fpath))
printl("""\nYou can edit this file before running 'barapost.py' in order to
  modify list of sequences that will be downloaded from Genbank
  and used as local (i.e. on your local computer) database by 'barapost.py'.""")
end_time = time()
printl('\n'+get_work_time() + " ({}) ".format(strftime("%Y-%m-%d %H:%M:%S", localtime(end_time))) + "- Probing task is completed\n")
logfile.close()
platf_depend_exit(0)