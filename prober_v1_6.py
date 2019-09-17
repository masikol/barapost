#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Version 1.6
# 08.09.2019 edition

# |===== Check python interpreter version =====|

from sys import version_info as verinf

if verinf.major < 3:#{
    print( "\nYour python interpreter version is " + "%d.%d" % (verinf.major, verinf.minor) )
    print("   Please, use Python 3!\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    raw_input("Press ENTER to exit:")
    exit(1)
#}

print("\n |=== prober.py (version 1.6) ===|\n")

# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, sleep, gmtime
start_time = time()


def get_work_time():#{
    return strftime("%H:%M:%S", gmtime( time() - start_time))
#}

# |===========================================|

import os
from re import search as re_search
from gzip import open as open_as_gzip # input files might be gzipped
from xml.etree import ElementTree # for retrieving information from XML BLAST report
from sys import intern

import http.client
import urllib.request
from urllib.error import HTTPError
import urllib.parse
import socket


# |===== Function that asks to press ENTER on Windows =====|

from sys import platform

def platf_depend_exit(exit_code):#{
    """
    Function asks to press ENTER press on Windows
        and exits after that.

    :type exit_code: int;
    """
    if platform.startswith("win"):
        input("Press ENTER to exit:")
    exit(exit_code)
#}


def print_error(text):#{
    """Function for printing error messages"""
    print("\n   \a!! - ERROR: " + text + '\n')
#}


from sys import stdout as sys_stdout
def printn(text):#{
    """
    Function prints text to the console without adding '\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
    """
    sys_stdout.write(text)
#}


# |===== Handle command line arguments =====|
help_msg = """
DESCRIPTION:\n
 prober.py -- script is designed for determinating the taxonomic position
    of nucleotide sequences by aligning each of them and regarding the best hit.\n
 The main goal of this script is to send a probing batch of sequences to BLAST server
    and discover, what Genbank records can be downloaded and used for further processing
    on your local machine by 'barapost.py'.\n
 This script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.\n
 Results of the work of this script are written to TSV files, that can be found in result directory:\n
  1. There is a file named "...probe_acc_list.tsv". It contains accessions and names of Genbank records that
    can be used for further processing on your local machine by 'barapost.py'.\n
  2. There is a file named "...result.tsv". It contains full result of aligning.
    Results of barapost.py's work will be appended to this file.\n
 FASTQ files processed by this script are meant to be processed afterwards by 'barapost.py'.
----------------------------------------------------------

Default parameters:\n
- all FASTQ and FASTA files in current directory will be processed;
- packet size (see '-p' option): 100 sequences;
- probing batch size (see '-b' option): 200 sequences;
- algorithm (see '-a' option): 'megaBlast';
- organisms (see '-g' option): full 'nt' database, i.e. no slices;
- output directory ('-o' option): directory named "barapost_result"
  nested in current directory;

  Default behavior of this script is to send certain batch (see '-b' option) of sequences to BLAST server.
It means that you should not process all your data by 'prober.py' -- it would take long time.\n
  Instead of this you should process some sequences by 'prober.py' -- it will determine,
what Genbank records (genomes, if you want) are present in your data and then go to 'barapost.py'.\n
  'barapost.py' will process the rest of you sequences in the same way like 'prober.py', but on your local computer.
'barapost.py' uses 'blast+' toolkit for this purpose. It would be much faster.\n
  Obviously, a probing batch cannot cover all variety of a data set,
so some sequences can be recognized as "unknown" while processing by 'barapost.py'.
But you always can run 'prober.py' again on "unknown" sequences.
----------------------------------------------------------

OPTIONS:\n
    -h (--help) --- show help message;\n
    -f (--infile) --- input FASTQ or FASTA (or '.fastq.gz', '.fasta.gz') file;
        You can specify multiple input files with this option (see EXAMPLES #2);\n
    -d (--indir) --- directory which contains FASTQ of FASTA files meant to be processed.
        I.e. all FASTQ and FASTA files in this direcory will be processed;
        Input files can be gzipped.\n
    -o (--outdir) --- output directory;\n
    -p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
        Value: integer number [1, 500]. Default value is 100;\n
    -a (--algorithm) --- BLASTn algorithm to use for aligning.
        Available values: 'megaBlast', 'discoMegablast', 'blastn'.
        Default is megaBlast;\n
    -g (--organisms) --- 'nt' database slices, i.e. organisms that you expect to see in result files.
        More clearly, functionality of this option is totally equal to "Organism" text boxes
        on this BLASTn page:
         'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'.
        Format of value: 
          <organism1_name>,<organism1_taxid>+<organism2_name>,<organism2_taxid>+...
        See EXAMPLES #2 and #3 below.
        Spaces are not allowed. Number of organisms can be from 1 to 5 (5 is maximum).
        Default value is full 'nt' database.
        You can find your Taxonomy IDs here: 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi'.\n
    -b (--probing-batch-size) --- number of sequences that will be aligned on BLAST server.
        After that a local database will be builded according to results of probing aligning.
        More clearly: records-hits will be downloaded from Genbank and will be used
        as local database. Further aligning against this database will be preformed
        on local machine with 'blast+' toolkit.
        Value: positive integer number. Default value is 200;\n
----------------------------------------------------------

EXAMPLES:\n
  1) Process one FASTQ file with default settings:\n
       ./prober.py -f reads.fastq\n
  2) Process FASTQ file and FASTA file with discoMegablast, packet size of 5 sequences.
     Search only among Erwinia sequences:\n
       ./prober.py -f reads.fastq.gz -f another_sequences.fasta -a discoMegablast -p 5 -g Erwinia,551\n
  3) Process all FASTQ and FASTA files in directory named "some_dir".
     Search only among Escherichia and viral sequences:\n
       ./prober.py -d some_dir -g Escherichia,561+viruses,10239 -o outdir
"""
from sys import argv
import getopt

try:#{
    opts, args = getopt.getopt(argv[1:], "hf:d:o:p:a:g:b:",
        ["help", "infile=", "indir=", "outdir=", "packet-size=", "algorithm=", "organisms=", "probing-batch-size="])
#}
except getopt.GetoptError as gerr:#{
    print( str(gerr) )
    platf_depend_exit(2)
#}

is_fq_or_fa = lambda f: True if not re_search(r".*\.(m)?f(ast)?(a|q)(\.gz)?$", f) is None else False

# Default values:
fq_fa_list = list()
indir_path = None
outdir_path = "prober_result"
packet_size = 100
probing_batch_size = 200
blast_algorithm = "megaBlast"
organisms = list() # default is whole 'nt' database

if len(args) != 0:#{
    print_error("prober.py does not take any positional arguments")
    print("Here are positional arguments specified by you:")
    for a in args:
        print(' ' + a)
    platf_depend_exit(1)
#}

for opt, arg in opts:#{
    
    if opt in ("-h", "--help"):#{
        print(help_msg)
        platf_depend_exit(0)
    #}

    if opt in ("-f", "--infile"):#{
        if not os.path.exists(arg):#{
            print_error("file '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        #}
        if not is_fq_or_fa(arg):#{
            print_error("file '{}' is not '.fastq' or '.fasta'!".format(arg))
            platf_depend_exit(1)
        #}
        fq_fa_list.append( os.path.abspath(arg) )
    #}

    if opt in ("-d", "--indir"):#{
        if not os.path.exists(arg):#{
            print_error("directory '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        #}
        if not os.path.isdir(arg):#{
            print_error("'{}' is not a directory!".format(arg))
            platf_depend_exit(1)
        #}
        indir_path = os.path.abspath(arg)

        fq_fa_list.extend(list( filter(is_fq_or_fa, os.listdir(indir_path)) ))
    #}

    if opt in ("-o", "--outdir"):#{
        outdir_path = os.path.abspath(arg)
    #}

    if opt in ("-p", "--packet-size"):#{
        try:#{
            packet_size = int(arg)
            if packet_size < 1 or packet_size > 500:
                raise ValueError
        #}
        except ValueError:#{
            print_error("packet_size (-p option) must be integer number from 1 to 500")
            platf_depend_exit(1)
        #}
    #}

    if opt in ("-a", "--algorithm"):#{
        if not arg in ("megaBlast", "discoMegablast", "blastn"):#{
            print_error("invalid value specified by '-a' option!")
            print("Available values: 'megaBlast', 'discoMegablast', 'blastn'")
            platf_depend_exit(1)
        #}
        blast_algorithm = arg
    #}

    if opt in ("-g", "--organisms"):#{
        max_org = 5

        try:#{
            org_list = arg.strip().split('+')
            org_list = list( map(str.strip, org_list) )

            if len(org_list) > max_org:#{
                raise Exception("\nYou can specify from 1 to {} organisms.\a".format(max_org))
            #}

            for org in org_list:#{
                name_and_taxid = org.strip().split(',')
                name_and_taxid = list( map(str.strip, name_and_taxid) )
                if len(name_and_taxid) != 2:#{
                    raise Exception("""\nOrganism's name and it's taxid should be separated by comma (,),
    and different organisms -- by plus (+).\n  Type for help: ./baparost.py -h\a""")
                #}
                # Validate TaxID integer format: it will raise ValueError if taxid is invalid
                tmp_taxid = int(name_and_taxid[1])
                if tmp_taxid < 1:#{
                    raise ValueError("\nTaxID should be positive integer number\a")
                #}
                organisms.append("{} (taxid:{})".format(name_and_taxid[0], name_and_taxid[1]))
            #}
        #}
        except ValueError:#{
            print("\n" + "=/"*20 + "\n")
            print_error("TaxID should be positive integer number\a")
            platf_depend_exit(1)
        #}
        except Exception as err:#{
            print("\n" + "=/"*20 + "\n")
            print_error("ERROR: invalid organisms (-g option) input format")
            print( str(err) )
            platf_depend_exit(1)
        #}
    #}

    if opt in ("-b", "--probing-batch-size"):#{
        try:#{
            probing_batch_size = int(arg)
            if probing_batch_size <= 0:
                raise ValueError
        #}
        except ValueError:#{
            print_error("probing batch size ('-b' option) must be positive integer number!")
            platf_depend_exit(1)
        #}
    #}
#}

# If no FASTQ or FASTA file have been found
if len(fq_fa_list) == 0:#{
    # If input directory was specified -- exit
    if not indir_path is None:#{
        print_error("""no input FASTQ or FASTA files specified
    or there is no FASTQ and FASTA files in the input directory.\n""")
        platf_depend_exit(1)
    #}
    # If input directory was not specified -- look for FASTQ files in current directory
    else:#{
        fq_fa_list = list( filter(is_fq_or_fa, os.listdir('.')) )
        if len(fq_fa_list) == 0:#{
            print_error("there is no FASTQ or FASTA files to process found.")
            platf_depend_exit(1)
        #}
    #}
#}

del help_msg # we do not need it any more

print( get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(start_time))) + "- Start working\n")


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

# |=== Delimiter for result tsv file ===|
DELIM = '\t'


# |=== File format constants ===|
FASTQ_LINES_PER_READ = 4
FASTA_LINES_PER_SEQ = 2

# |=== Check if there is enough sequeneces in files (>= probing_batch_size) ===|
seqs_at_all = 0
print() # just print new line
for file in fq_fa_list:#{
    how_to_open = OPEN_FUNCS[ is_gzipped(file) ]
    if is_fastq(file):
        n_seqs = int( sum(1 for line in how_to_open(file, 'r')) / FASTQ_LINES_PER_READ )
    else:
        n_seqs = sum(1 if line[0] == '>' else 0 for line in how_to_open(file, 'r'))
    seqs_at_all += n_seqs
    if seqs_at_all >= probing_batch_size:
        break
#}

# Print a warning message if a user has specified batch size that is greater than number of sequences he has at all
# And do not disturb him if he ran 'prober.py' with default batch size
if seqs_at_all < probing_batch_size and ("-b" in argv or "--probing_batch_size" in argv):#{
    print('\n'+'-'*20)
    print("\a  Warning!\n There are totally {} sequences in your files.".format(seqs_at_all))
    print(" Probing batch size specified by you is {}".format(probing_batch_size))

    while True:#{
        reply = input("\nPress ENTER to process all your sequences anyway.\n  Or enter 'q' to exit:>>")
        if reply == "":
            probing_batch_size = seqs_at_all
            break
        elif reply == 'q':
            platf_depend_exit(0)
        else:
            print_error("invalid reply")
            continue
    #}
#}

if seqs_at_all < probing_batch_size and not ("-b" in argv or "--probing_batch_size" in argv):#{
    probing_batch_size = seqs_at_all
#}


# |===== Function for checking if 'https://blast.ncbi.nlm.nih.gov' is available =====|

def check_connection():#{
    """
    Function checks if 'https://blast.ncbi.nlm.nih.gov' is available.

    :return: None if 'https://blast.ncbi.nlm.nih.gov' is available;
    """

    try:#{

        ncbi_server = "https://blast.ncbi.nlm.nih.gov"
        status_code = urllib.request.urlopen(ncbi_server).getcode()

        # Just in case
        if status_code != 200:#{
            print('\n' + get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            print("Check your Internet connection.\a")
            print("Status code: {}".format(status_code))
            platf_depend_exit(-2)
        #}
        return
    #}
    except OSError as err:#{

        print('\n' + get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
        print("Check your Internet connection.\a")
        print('\n' + '=/' * 20)
        print( str(err) )

        # 'urllib.request.HTTPError' can provide a user with information about the error
        if isinstance(err, HTTPError):#{
            print("Status code: {}".format(err.code))
            print(err.reason)
        #}
        platf_depend_exit(-2)
    #}
#}


# |===== Question funtions =====|

def is_continued():#{
    """
    Function asks the user if he/she wants to continue the previous run.

    :return: True if the decision is to continue, else False;
    :return type: bool;
    """

    continuation = None

    while continuation is None:#{
        continuation = input("""
Would you like to continue the previous run?
    1. Continue!
    2. Start from the beginning.

Enter the number (1 or 2):>> """)
        # Check if entered value is integer number. If no, give another attempt.
        try:#{
            continuation = int(continuation)
            if continuation != 1 and continuation != 2:#{ Check if input number is 1 or 2
                print("\n   Not a VALID number entered!\a\n" + '~'*20)
                continuation = None
            #}
            else:#{
                print("You have chosen number " + str(continuation) + '\n')
            #}
        #}
        except ValueError:#{
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            continuation = None
        #}
    #}
    return(True if continuation == 1 else False)
#}


def get_packet_size(num_reads):#{
    """
    Function asks the user about how many query sequences will be sent 
        to NCBI BLAST as a particular request.

    :return: the number of query sequences;
    :return type: int;
    """

    packet_size = None
    # You cannot sent more query sequences than you have
    limit = num_reads if num_reads <= 500 else 500

    while packet_size is None:#{
        
        packet_size = input("""
Please, specify the number of sequences that should be sent to the NCBI server in one request.
E.g. if you have 10 sequences in your file, you can send 10 sequences as single
    request -- in this case you should enter number 10. You may send 2 requests containing
    5 sequences both -- in this case you should enter number 5.


There are {} sequences left to process in current file.
Enter the number (from 1 to {}):>> """.format(num_reads, limit))
        # Check if entered value is integer number. If no, give another attempt.
        try:#{
            packet_size = int(packet_size)
            if packet_size < 1 or packet_size > limit:#{ Check if input number is in [1, limit]
                print("\n   Not a VALID number entered!\a\n" + '~'*20)
                packet_size = None
            #}
            else:#{
                print("You have chosen number " + str(packet_size) + '\n')
            #}
        #}
        except ValueError:#{
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            packet_size = None
        #}
    #}
    return(packet_size)
#}

# |===== End of question funtions =====|


get_phred33 = lambda q_symb: ord(q_symb) - 33

def get_read_avg_qual(qual_str):#{

    phred33 = map(get_phred33, list(qual_str))
    read_qual = round( sum(phred33) / len(qual_str), 2 )
    return read_qual
#}

def configure_qual_dict(fastq_path):#{

    qual_dict = dict()
    how_to_open = OPEN_FUNCS[ is_gzipped(fastq_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq_path) ]

    with how_to_open(fastq_path) as fastq_file:#{
        counter = 1
        line = fmt_func(fastq_file.readline())
        while line != "":#{
            if counter == 1:#{
                seq_id = intern( (line.partition(' ')[0])[1:] )
            #}
            counter += 1
            line = fmt_func(fastq_file.readline())
            if counter == 4:#{
                qual_dict[seq_id] = get_read_avg_qual(line)
                counter = 0
            #}
        #}
    #}

    return qual_dict
#}


def fastq2fasta(fq_fa_path, i, new_dpath):#{
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
    
    fasta_path = os.path.basename(fq_fa_path).replace(".fastq", ".fasta") # change extention
    fasta_path = fasta_path.replace(".gz", "") # get rid of ".gz" extention
    fasta_path = os.path.join(new_dpath, fasta_path) # place FASTA file into result directory

    how_to_open = OPEN_FUNCS[ is_gzipped(fq_fa_path) ]

    fastq_patt = r".*\.f(ast)?q(\.gz)?$"

    num_lines = 0 # variable for counting lines in a file
    if not re_search(fastq_patt, fq_fa_path) is None:#{

        global FASTQ_LINES_PER_READ

        # Get ready to process gzipped files
        # how_to_open = OPEN_FUNCS[ is_gzipped(fq_fa_path) ]
        fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_fa_path) ]

        with how_to_open(fq_fa_path) as fastq_file, open(fasta_path, 'w') as fasta_file:#{

            counter = 1 # variable for retrieving only 1-st and 2-nd line of FASTQ record
            for line in fastq_file:#{
                line = fmt_func(line)
                if counter <= 2:#{      write only 1-st and 2-nd line out of 4
                    if line[0] == '@':
                        line = '>' + line[1:]  # replace '@' with '>'
                    fasta_file.write(line + '\n')
                #}
                # reset the counter if the 4-th (quality) line has been encountered
                elif counter == 4:
                    counter = 0
                counter += 1
                num_lines += 1
            #}
        #}
        num_reads = int(num_lines / FASTQ_LINES_PER_READ) # get number of sequences

        print("\n{}. '{}' ({} reads) --> FASTA".format(i+1, os.path.basename(fq_fa_path), num_reads))
    #}
    # IF FASTA file is already created
    # We need only number of sequences in it.
    elif not re_search(fastq_patt, fq_fa_path) is None and os.path.exists(fasta_path):#{
        num_lines = sum(1 for line in how_to_open(fq_fa_path)) # get number of lines
        num_reads = int( num_lines / FASTQ_LINES_PER_READ ) # get number of sequences
    #}
    # We've got FASTA source file
    # We need only number of sequences in it.
    else:#{ 
        global FASTA_LINES_PER_SEQ

        # num_lines = sum(1 for line in how_to_open(fq_fa_path)) # get number of lines
        # num_reads = int( num_lines / FASTA_LINES_PER_SEQ ) # get number of sequences
        num_reads = sum(1 if line[0] == '>' else 0 for line in how_to_open(fq_fa_path, 'r'))
        fasta_path = fq_fa_path
    #}

    print("\n |===== file: '{}' ({} sequences) =====|".format(os.path.basename(fasta_path), num_reads))
    return {"fpath": fasta_path, "nreads": num_reads}
#}


def rename_file_verbosely(file, directory):#{

    if os.path.exists(file):#{

        is_analog = lambda f: file[:os.path.basename(file).rfind('.')] in f
        num_analog_files = len( list(filter(is_analog, os.listdir(directory))) )

        print('\n' + get_work_time() + " - Renaming old file:")
        name_itself = file[: file.rfind('.')]
        ext = file[file.rfind('.'):]
        num_analog_files = str(num_analog_files)
        new_name = name_itself+"_old_"+num_analog_files+ext

        print("   '{}' --> '{}'".format(os.path.basename(file), new_name))
        os.rename(file, new_name)
    #}
#}


def look_around(outdir_path, new_dpath, fasta_path, blast_algorithm):#{
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
    tsv_res_fpath = "{}_{}_result.tsv".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)
    # Form path to accession file
    acc_fpath = os.path.join(outdir_path, "{}_probe_acc_list.tsv".format(blast_algorithm))

    num_done_reads = None # variable to keep number of succeffdully processed sequences

    continuation = None
    # Check if there are results from previous run.
    if os.path.exists(tsv_res_fpath) or os.path.exists(tmp_fpath):#{
        print('\n' + get_work_time() + " - The previous result file is found in the directory:")
        print("   '{}'".format(new_dpath))
        # Allow politely to continue from last successfully sent packet.
        continuation = is_continued()
        if not continuation:
            rename_file_verbosely(tsv_res_fpath, new_dpath)
            rename_file_verbosely(tmp_fpath, new_dpath)
            rename_file_verbosely(acc_fpath, new_dpath)
    #}

    if continuation == True:#{   Find the name of last successfull processed sequence
        print("Let's try to continue...")

        # Collect information from result file
        if os.path.exists(tsv_res_fpath):#{
            try:#{ There can be invalid information in this file
                with open(tsv_res_fpath, 'r') as res_file:#{
                    lines = res_file.readlines()
                    num_done_reads = len(lines) - 1 # the first line is a head
                    last_line = lines[-1]
                    last_seq_id = last_line.split(DELIM)[0]
                #}
            except Exception as err:#{
                print("\nData in result file '{}' not found or broken.".format(tsv_res_fpath))
                print( str(err) )
                print("Start from the beginning.")
                rename_file_verbosely(tsv_res_fpath, new_dpath)
                rename_file_verbosely(tmp_fpath, new_dpath)
                rename_file_verbosely(acc_fpath, new_dpath)
                return None
            #}
            else:#{
                print("Last sent sequence: " + last_seq_id)
        #}


        # Collect information from accession file
        global acc_dict
        if os.path.exists(acc_fpath):#{
            try:#{  There can be invalid information in this file
                with open(acc_fpath, 'r') as acc_file:#{
                    lines = acc_file.readlines()[9:] # omit description and head of the table
                    local_files_filtered = list( filter(lambda x: False if os.path.exists(x) else True, lines) ) # omit file paths
                    for line in local_files_filtered:#{
                        vals = line.split(DELIM)
                        acc = intern(vals[0].strip())
                        acc_dict[acc] = (vals[1].strip(), vals[2].strip(), vals[3].strip())
                    #}
                #}
            #}
            except Exception as err:#{
                print("\nData in accession file '{}' not found or broken.".format(acc_fpath))
                print( str(err) )
                print("Start from the beginning.")
                rename_file_verbosely(tsv_res_fpath, new_dpath)
                rename_file_verbosely(tmp_fpath, new_dpath)
                rename_file_verbosely(acc_fpath, new_dpath)
                return None
            #}
            else:#{
                print("\nHere are Genbank records encountered during previous run:")
                for acc in acc_dict.keys():#{
                    s_letter = "s" if acc_dict[acc][2] > 1 else ""
                    print(" {} sequence{} - {}, '{}'".format(acc_dict[acc][2], s_letter, acc, acc_dict[acc][1]))
                #}
                print()
            #}
        #}

        # If we start from the beginning, we have no sequences processed
        if num_done_reads is None:
            num_done_reads = 0

        # Get packet size, number of the last sent packet and RID
        try:#{ There can be invalid information in tmp file of tmp file may not exist
            with open(tmp_fpath, 'r') as tmp_file:
                temp_lines = tmp_file.readlines()
            packet_size = int(re_search(r"packet_size: ([0-9]+)", temp_lines[0]).group(1).strip())
            saved_npack = int(re_search(r"sent_packet_num: ([0-9]+)", temp_lines[1]).group(1).strip())
            RID_save = re_search(r"Request_ID: (.+)", temp_lines[2]).group(1).strip()
            # If aligning is performed on local machine, there is no reason for requesting results.
            # Therefore this packet will be aligned once again.
        #}
        except Exception as exc:#{
            total_num_seqs = sum(1 if line[0] == '>' else 0 for line in how_to_open(fasta_path))
            print("\n    Attention! Temporary file not found or broken!")
            print( str(exc) )
            print("{} sequences have been already processed.".format(num_done_reads))
            print("There are {} sequences totally in file '{}'.".format(total_num_seqs, os.path.basename(fasta_path)))
            print("Maybe you've already processed this file with 'prober.py'?\n")

            reply = "BULLSHIT"
            global omit_file
            while True:#{
                reply = input("Press ENTER to omit this file and go to the next one.\n \
 Or enter 'c' to to continue processing this file:>>")
                if reply == 'c':
                    if num_done_reads >= total_num_seqs:#{
                        print_error("there are no sequences left to process in this file")
                        print("Omitting this file")
                        omit_file = True
                        return
                    #}
                    break
                elif reply == "":
                    omit_file = True
                    global seqs_processed
                    seqs_processed += num_done_reads
                    return
                else:
                    print_error("invalid reply")
            #}


            print("{} reads have been already processed".format(num_done_reads))
            print("{} reads left".format(total_num_seqs - num_done_reads))
            packet_size = get_packet_size(total_num_seqs - num_done_reads)
            return {
                "pack_size": packet_size,
                "sv_npck": int(num_done_reads / packet_size),
                "RID": None,
                "acc_fpath": acc_fpath,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        else:#{
            # Return data from previous run
            return {
                "pack_size": packet_size,
                "sv_npck": saved_npack,
                "RID": RID_save,
                "acc_fpath": acc_fpath,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        #}
    elif continuation == False:#{
        # If we've decided to start from the beginnning - rename old files
        rename_file_verbosely(tsv_res_fpath, new_dpath)
        rename_file_verbosely(tmp_fpath, new_dpath)
        rename_file_verbosely(acc_fpath, new_dpath)
    #}
    return None
#}


def get_packet(fasta_file, packet_size, fmt_func):#{
    """
    Function collects the packet of query sequences to send to BLAST server.

    :param fasta_file: file instance of FASTA file to retrieve sequences from;
    :type fasta_file: _io.TextIOWrapper;
    :param packet_size: number of query sequences to send as a particular request;
    :type packet_size: int;

    Returns dict of the following structure:
    {
        "fasta": FASTA_data_containing_query_sequences (str),
        "names": list_of_sequence_ids_from_fasta_file (list<str>)
    }
    """

    global next_id_line
    line = fmt_func(fasta_file.readline())
    if line.startswith('>') and ' ' in line:#{:
        line = line.partition(' ')[0]+'\n'
    #}
    packet = ""

    i = 0
    if not next_id_line is None:
        packet += next_id_line
    packet += line

    while i < packet_size:#{

        line = fmt_func(fasta_file.readline())
        if line.startswith('>'):#{
            if ' ' in line:
                line = line.partition(' ')[0]+'\n'
            i += 1
        #}
        if line == "":
            break
        packet += line
    #}

    if line != "":#{
        next_id_line = packet.splitlines()[-1]+'\n'
        packet = '\n'.join(packet.splitlines()[:-1])
    #}
    else:#{
        next_id_line = None
        packet = packet.strip()
    #}

    names = list( filter(lambda l: True if l.startswith('>') else False, packet.splitlines()) )
    names = list( map(lambda l: l.partition(' ')[0].strip(), names) )

    return {"fasta": packet.strip(), "names": names}
#}


def remove_tmp_files(*paths):#{
    """
    Function removes files passed to it.
    Actually, passed arguments are paths ('str') to files meant to be removed.
    """
    for path in paths:#{
        if os.path.exists(path):#{
            os.unlink(path)
        #}
    #}
#}


def configure_request(packet, blast_algorithm, organisms):#{
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

    # 'nt' database slices:
    for i, org in enumerate(organisms):#{
        payload["EQ_MENU{}".format(i if i > 0 else "")] = org
    #}
    payload["NUM_ORG"] = str( len(organisms) )

    payload = urllib.parse.urlencode(payload)
    headers = { "Content-Type" : "application/x-www-form-urlencoded" }

    return {"payload":payload, "headers": headers}
#}

    
def send_request(request, pack_to_send, packs_at_all, filename, tmp_fpath):#{
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

    # Save packet size
    with open(tmp_fpath, 'w') as tmp_file:
        tmp_file.write("packet_size: {}\n".format(packet_size))

    while error:#{
        try:#{
            conn = http.client.HTTPSConnection(server) # create a connection
            conn.request("POST", url, payload, headers) # send the request
            response = conn.getresponse() # get the response
            response_text = str(response.read(), "utf-8") # get response text
        #}
        except OSError as oserr:#{
            print(get_work_time() + "\n - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            print( repr(err) )
            print("barapost will try to connect again in 30 seconds...\n")
            sleep(30)
        #}
        else:#{ if no exception occured
            error = False
        #}
    #}

    try:#{
        rid = re_search(r"RID = (.+)", response_text).group(1) # get Request ID
        rtoe = int(re_search(r"RTOE = ([0-9]+)", response_text).group(1)) # get time to wait provided by the NCBI server
    #}
    except AttributeError:#{
        print_error("seems, ncbi has denied your request.")
        print("Response is in file 'request_denial_response.html'")
        with open("request_denial_response.html", 'w') as den_file:#{
            den_file.write(response_text)
        #}
        exit(1)
    #}
    finally:#{
        conn.close()
    #}

    # Save temporary data
    with open(tmp_fpath, 'a') as tmpfile:
        tmpfile.write("sent_packet_num: {}\n".format(pack_to_send))
        tmpfile.write("Request_ID: {}".format(rid))

    # /=== Wait for alignment results ===\

    return( wait_for_align(rid, rtoe, pack_to_send, packs_at_all, filename) )
#}


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

    print("\n{} - Requesting for alignment results: {}, '{}' ({}/{})".format(get_work_time(),
    rid, filename, pack_to_send, packs_at_all))
    if rtoe > 0:#{ RTOE can be zero at the very beginning of continuation
        print("{} - BLAST server estimates that alignment will be accomplished in {} seconds ".format(get_work_time(), rtoe))
        print("{} - Waiting for {}+3 (+3 extra) seconds...".format(get_work_time(), rtoe))
        # Server migth be wrong -- we will give it 3 extra seconds
        sleep(rtoe + 3)
        print("{} - {} seconds have passed. Checking if alignment is accomplished...".format(get_work_time(), rtoe+3))
    #}

    server = "blast.ncbi.nlm.nih.gov"
    wait_url = "/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid
    there_are_hits = False

    while True:#{
        error = True
        while error:#{
            try:
                conn = http.client.HTTPSConnection(server) # create connection
                conn.request("GET", wait_url) # ask for if there areresults
            except TimeoutError as err:
                print("{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
                print("   " + str(err))
                error = True
                sleep(30)
            except http.client.RemoteDisconnected as err:
                print("{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
                print("   " + str(err))
                error = True
                sleep(30)
            except socket.gaierror as err:
                print("{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
                print("   " + str(err))
                error = True
                sleep(30)
            except ConnectionResetError as err:
                print("{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
                print("   " + str(err))
                error = True
                sleep(30)
            except FileNotFoundError as err:
                print("{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
                print("   " + str(err))
                error = True
                sleep(30)
            except http.client.CannotSendRequest as err:
                print("{} - Unable to connect to the NCBI server. Let's try to connect in 30 seconds.".format(get_work_time()))
                print("   " + str(err))
                error = True
                sleep(30)
            else:
                error = False # if no exception ocured
        #}

        response = conn.getresponse() # get the resonse
        resp_content = str(response.read(), "utf-8") # get response text
        conn.close()

        if not re_search("Status=WAITING", resp_content) is None:#{ if server asks to wait
            printn("\r{} - The request is still processing. Waiting for 60 seconds".format(get_work_time()))
            for i in range(1, 7):#{ indicate each 10 seconds with a dot
                sleep(10)
                printn("\r{} - The request is still processing. Waiting for 60 seconds{}".format(get_work_time(), '.'*i))
            #}
            print() # go to next line
            continue
        #}
        if not re_search("Status=FAILED", resp_content) is None:#{ if job failed
            print('\n' + get_work_time() + " - Job failed\a\n")
            response_text = """{} - Job for query {} ({}/{}) with Request ID {} failed.
    Contact NCBI or try to start it again.\n""".format(get_work_time(), filename, pack_to_send, packs_at_all, rid)
            return None
        #}
        if not re_search("Status=UNKNOWN", resp_content) is None:#{ if job is expired
            print('\n' + get_work_time() + " - Job expired\a\n")
            respond_text = """{} - Job for query {} ({}/{}) with Request ID {} expired.
    Try to start it again\n""".format(get_work_time(), filename, pack_to_send, packs_at_all, rid)
            return "expired"
        #}
        if not re_search("Status=READY", resp_content) is None:#{ if results are ready
            there_are_hits = True
            print("\n{} - Result for query '{}' ({}/{}) is ready!".format(get_work_time(), filename, pack_to_send, packs_at_all))
            if re_search("ThereAreHits=yes", resp_content) is not None:#{ if there are hits
                for i in range(45, 0, -5):
                    print('-' * i)
                print() # just print a blank line
                break
            #}
            else:#{ if there are no hits
                print(get_work_time() + " - There are no hits. It happens.\n")
                break
            #}
        #}
        # Execution should not reach here
        print('\n' + get_work_time() + " - Fatal error. Please contact the developer.\a\n")
        platf_depend_exit(1)
    #}

    # Retrieve XML result
    retrieve_xml_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&ALIGNMENTS=1&RID=" + rid
    conn = http.client.HTTPSConnection(server)
    conn.request("GET", retrieve_xml_url)
    response = conn.getresponse()

    respond_text = str(response.read(), "utf-8")
    conn.close()

    # Retrieve human-readable text and put it into result directory
    if there_are_hits:#{
        save_txt_align_result(server, filename, pack_to_send, rid)
    #}

    return respond_text
#}


def save_txt_align_result(server, filename, pack_to_send, rid):#{

    global outdir_path

    retrieve_text_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&DESCRIPTIONS=1&ALIGNMENTS=1&RID=" + rid
    conn = http.client.HTTPSConnection(server)
    conn.request("GET", retrieve_text_url)
    response = conn.getresponse()

    txt_hpath = os.path.join(outdir_path, filename.replace(".fasta", ""), "blast_result_{}.txt".format(pack_to_send))
    # Write text result for a human to read
    with open(txt_hpath, 'w') as txt_file:
        txt_file.write(str(response.read(), "utf-8") + '\n')
    conn.close()
#}


def parse_align_results_xml(xml_text, seq_names):#{
    """
    Function parses BLAST xml response and returns tsv lines containing gathered information:
        1. Query name.
        2. Hit name formatted by 'format_taxonomy_name()' function.
        3. Hit accession.
        4. Length of alignment.
        5. Percent of identity.
        6. Percent of gaps.
        7. E-value.
    Erroneous tsv lines that function may produce:
        1. "<query_name>\\tQuery has been lost: ERROR, Bad Gateway"
            if data packet has been lost.
        2. "<query_name>\\tQuery has been lost: BLAST ERROR"
            if BLAST error occured.
        3. "<query_name>\\tNo significant similarity found"
            if no significant similarity has been found
        Type of return object: list<str>.
    """

    result_tsv_lines = list()

    # /=== Validation ===/

    if "Bad Gateway" in xml_text:#{
        print('\n' + '=' * 45)
        print(get_work_time() + " - ERROR! Bad Gateway! Data from last packet has lost.")
        print("It would be better if you restart the script.")
        print("Here are names of lost queries:")
        for i, name in enumerate(seq_names):#{
            print("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + DELIM + "Query has been lost: ERROR, Bad Gateway")
        #}
        input("Press ENTER to continue...")

        return result_tsv_lines
    #}

    if "to start it again" in xml_text:#{
        print('\n' + get_work_time() + "BLAST ERROR!")

        print("Here are names of lost queries:")
        for i, name in enumerate(seq_names):#{
            print("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + DELIM +"Query has been lost: BLAST ERROR")
        #}

        input("Press ENTER to continue...")

        return result_tsv_lines
    #}

    # /=== Parse BLAST XML response ===/

    root = ElementTree.fromstring(xml_text) # get tree instance

    global new_acc_dict
    global qual_dict

    # Iterate through "Iteration" and "Iteration_hits" nodes
    for iter_elem, iter_hit in zip(root.iter("Iteration"), root.iter("Iteration_hits")):#{
    
        # "Iteration" node contains query name information
        query_name = intern(iter_elem.find("Iteration_query-def").text)

        query_len = iter_elem.find("Iteration_query-len").text

        if not qual_dict is None:#{
            ph33_qual = qual_dict[query_name]
            miscall_prop = round(10**(ph33_qual/-10), 3)
            accuracy = round( 100*(1 - miscall_prop), 2 ) # expected percent of correctly called bases
            qual_info_to_print = "    Average quality of this read is {}, i.e. accuracy is {}%;\n".format(ph33_qual,
                accuracy)
        #}
        else:#{
            # If FASTA file is processing, print dashed in quality columns
            ph33_qual = "-"
            accuracy = "-" # expected percent of correctly called bases
            qual_info_to_print = ""
        #}

        # If there are any hits, node "Iteration_hits" contains at least one "Hit" child
        hit = iter_hit.find("Hit")
        if hit is not None:#{

            # Get full hit name (e.g. "Erwinia amylovora strain S59/5, complete genome")
            hit_name = hit.find("Hit_def").text
            # Format hit name (get rid of stuff after comma)
            hit_taxa_name = hit_name[: hit_name.find(',')] if ',' in hit_name else hit_name
            hit_taxa_name = hit_taxa_name.replace(" complete genome", "") # sometimes there are no comma before it
            hit_taxa_name = hit_taxa_name.replace(' ', '_')


            hit_acc = intern(hit.find("Hit_accession").text) # get hit accession
            gi_patt = r"gi\|([0-9]+)" # pattern for GI number finding
            hit_gi = re_search(gi_patt, hit.find("Hit_id").text).group(1)
            try:#{
                new_acc_dict[hit_acc][2] += 1
            #}
            except KeyError:#{
                new_acc_dict[hit_acc] = [hit_gi, hit_name, 1]
            #}


            # Find the first HSP (we need only the first one)
            hsp = next(hit.find("Hit_hsps").iter("Hsp"))

            align_len = hsp.find("Hsp_align-len").text.strip()

            pident = hsp.find("Hsp_identity").text # get number of matched nucleotides

            gaps = hsp.find("Hsp_gaps").text # get number of gaps

            evalue = hsp.find("Hsp_evalue").text # get e-value

            pident_ratio = round( float(pident) / int(align_len) * 100, 2)
            gaps_ratio = round( float(gaps) / int(align_len) * 100, 2)

            print("""\n '{}' -- '{}';
    Query length - {} nt;
    Identity - {}/{} ({}%); Gaps - {}/{} ({}%);""".format(query_name, hit_taxa_name,
                    query_len, pident, align_len, pident_ratio, gaps, align_len, gaps_ratio))

            # Append new tsv line containing recently collected information
            result_tsv_lines.append( DELIM.join( (query_name, hit_taxa_name, hit_acc, query_len,
                align_len, pident, gaps, evalue, str(ph33_qual), str(accuracy)) ))
        #}
        else:#{
            # If there is no hit for current sequence
            print("\n '{}' -- No significant similarity found;\n    Query length - {};".format(query_name, query_len))
            result_tsv_lines.append(DELIM.join( (query_name, "No significant similarity found", "-", query_len,
                "-", "-", "-", "-", str(ph33_qual), str(accuracy)) ))
        #}
        printn(qual_info_to_print)
    #}
    return result_tsv_lines
#}


def write_result(res_tsv_lines, tsv_res_path, acc_file_path, fasta_hname, outdir_path):#{
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
    if not os.path.exists(tsv_res_path):#{
        with open(tsv_res_path, 'w') as tsv_res_file:#{
            tsv_res_file.write(DELIM.join( ["QUERY_ID", "HIT_NAME", "HIT_ACCESSION", "QUERY_LENGTH",
                "ALIGNMENET_LENGTH", "IDENTITY", "GAPS", "E-VALUE", "AVG_PHRED33", "ACCURACY(%)"] ) + '\n')
        #}
    #}
    # Write reslut tsv lines to this file
    with open(tsv_res_path, 'a') as tsv_res_file:#{
        for line in res_tsv_lines:
            tsv_res_file.write(line + '\n')
    #}

    # === Write accession information ===

    global acc_dict
    global new_acc_dict
    global blast_algorithm
    acc_file_path = os.path.join(outdir_path, "{}_probe_acc_list.tsv".format(blast_algorithm))

    with open(acc_file_path, 'w') as acc_file:#{
        acc_file.write("# Here are accessions, GI numbers and descriptions of Genbank records that can be used for sorting by 'barapost.py'\n")
        acc_file.write("# Values in this file are delimited by tabs.\n")
        acc_file.write("# You are welcome to edit this file by adding,\n")
        acc_file.write("#   removing or muting lines (with adding '#' symbol in it's beginning, just like this description).\n")
        acc_file.write("# Lines muted with '#' won't be noticed by 'barapost.py'.\n")
        acc_file.write("# You can specify your own FASTA files that you want to use as database for 'barapost.py'.\n")
        acc_file.write("# To do it, just write your FASTA file's path to this TSV file in new line.\n\n")
        acc_file.write(DELIM.join( ["ACCESSION", "GI_NUMBER", "RECORD_NAME", "OCCURRENCE_NUMBER"] ) + '\n')
    #}

    for acc, other_info in new_acc_dict.items():#{
        try:#{
            acc_dict[acc][2] += 1
        #}
        except KeyError:#{
            acc_dict[acc] = other_info
        #}
    #}
    # Write accessions and record names
    with open(acc_file_path, 'a') as acc_file:#{
        for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][2]):
            acc_file.write(DELIM.join( (acc, other_info[0], other_info[1], str(other_info[2]))) + '\n')
        #}
    #}

    new_acc_dict.clear()
#}


def create_result_directory(fq_fa_path, outdir_path):#{
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
    if not os.path.exists(new_dpath):#{
        os.makedirs(new_dpath)
    #}

    return new_dpath
#}



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
#        "pack_size": packet_size (int),
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

print("Checking Internet connection...")
check_connection()
print("OK\n\n")


print(" - Probing batch size: {} sequences;".format(probing_batch_size))
print(" - Packet size: {} sequences;".format(packet_size))
print(" - BLAST algorithm: {};".format(blast_algorithm))
print(" - Database: nt;")
if len(organisms) > 0:#{
    for db_slice in organisms:
        print("   {};".format(db_slice))
#}
print()

print(" Following files will be processed:")
for i, path in enumerate(fq_fa_list):#{
    print("    {}. '{}'".format(i+1, path))
#}
print('-'*30 + '\n')

# Variable for counting accessions of records menat to be downloaded from Genbank.
# Is used only for printing the list of accessions to console.
acc_counter = 0
# Dictionary of accessions and record names.
# Accessions are keys, record names are values.
# This dictionary is filled while processing and at the beginning of continuation.
acc_dict = dict()
# Dictionary of accessions and record names encountered while sending of the current packet.
# Accessions are keys, record names are values.
new_acc_dict = dict()

# Counter of processed sequences
seqs_processed = 0

# Variable that contains id of next sequence in current FASTA file.
# If no or all sequences in current FASTA file have been already processed, this variable is None
# Function 'get_packet' changes this variable
next_id_line = None

stop = False
omit_file= False


# Iterate through found source FASTQ and FASTA files
for i, fq_fa_path in enumerate(fq_fa_list):#{

    # Configure quality dictionary
    qual_dict = configure_qual_dict(fq_fa_path) if is_fastq(fq_fa_path) else None
    
    # Create the result directory with the name of FASTQ of FASTA file being processed:
    new_dpath = create_result_directory(fq_fa_path, outdir_path)

    # Convert FASTQ file to FASTA (if it is FASTQ) and get it's path and number of sequences in it:
    curr_fasta = fastq2fasta(fq_fa_path, i, new_dpath)

    # "hname" means human readable name (i.e. without file path and extention)
    fasta_hname = os.path.basename(curr_fasta["fpath"]) # get rid of absolure path
    fasta_hname = re_search(r"(.*)\.(m)?f(ast)?a", fasta_hname).group(1) # get rid of file extention

    # Look around and ckeck if there are results of previous runs of this script
    # If 'look_around' is None -- there is no data from previous run
    previous_data = look_around(outdir_path, new_dpath, curr_fasta["fpath"],
        blast_algorithm)

    if omit_file:#{
        omit_file = False
        continue
    #}

    if previous_data is None:#{ # If there is no data from previous run

        num_done_reads = 0 # number of successfully processed sequences
        saved_npack = None # number of last sent packet (there is no such stuff for de novo run)
        tsv_res_path = "{}_{}_result.tsv".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form result tsv file path
        tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form temporary file path
        acc_fpath = os.path.join(outdir_path, "{}_probe_acc_list.tsv".format(blast_algorithm)) # form path to accession file
    #}
    else:#{ # if there is data from previous run

        num_done_reads = previous_data["n_done_reads"] # get number of successfully processed sequences
        saved_npack = previous_data["sv_npck"] # get number of last sent packet
        packet_size = previous_data["pack_size"] # packet size sholud be the same as it was in previous run
        tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        tmp_fpath = previous_data["tmp_fpath"] # temporary file sholud be the same as during previous run
        acc_fpath = previous_data["acc_fpath"] # accession file sholud be the same as during previous run
        saved_RID = previous_data["RID"] # having this RID we can try to get response for last request
        contin_rtoe = 0 # we will not sleep at the very beginning of continuation
    #}

    if num_done_reads != 0:
        seqs_processed = num_done_reads

    if (probing_batch_size - seqs_processed) <= curr_fasta["nreads"]:
        tmp_num = (probing_batch_size - seqs_processed)
    else:
        tmp_num = curr_fasta["nreads"]
    packs_at_all = tmp_num // packet_size # Calculate total number of packets sent from current FASTA file

    if tmp_num % packet_size > 0: # And this is ceiling (in order not to import 'math')
        packs_at_all += 1
    packs_received = int( num_done_reads / packet_size ) # number of successfully processed sequences

    if is_gzipped(curr_fasta["fpath"]):#{
        fmt_func = lambda l: l.decode("utf-8")
        how_to_open = open_as_gzip
    #}
    else:#{
        fmt_func = lambda l: l
        how_to_open = open
    #}

    with how_to_open(curr_fasta["fpath"]) as fasta_file:#{

        # Go untill the last processed sequence
        # for _ in range( int(num_done_reads * 2) ):#{
        #     fasta_file.readline()
        # #}
        get_packet(fasta_file, num_done_reads, fmt_func)

        reads_left = curr_fasta["nreads"] - num_done_reads # number of sequences left to precess
        packs_left = packs_at_all - packs_received # number of packets left to send
        pack_to_send = packs_received+1 if packs_received > 0 else 1 # number of packet meant to be sent now

        # Iterate througth packets left to send
        for i in range(packs_left):#{

            packet = get_packet(fasta_file, packet_size, fmt_func) # form the packet

            if packet["fasta"] is "":#{   Just in case
                print("Recent packet is empty")
                break
            #}

            print("\nGo to BLAST (" + blast_algorithm + ")!")
            print("Request number {} out of {}.".format(pack_to_send, packs_at_all))

            send = True

            # If current packet has been already send, we can try to get it and not to send it again
            if pack_to_send == saved_npack and saved_RID is not None:#{

                align_xml_text = wait_for_align(saved_RID, contin_rtoe,
                    pack_to_send, packs_at_all, fasta_hname+".fasta") # get BLAST XML response

                # If request is not expired get he result and not send it again
                if align_xml_text != "expired":#{
                    send = False

                    result_tsv_lines = parse_align_results_xml(align_xml_text,
                        packet["names"]) # get result tsv lines

                    seqs_processed += len( packet["names"] )

                    # Write the result to tsv
                    write_result(result_tsv_lines, tsv_res_path, acc_fpath, fasta_hname, outdir_path)
                #}
            #}

            if send:#{
            
                request = configure_request(packet["fasta"], blast_algorithm, organisms) # get the request

                # Send the request get BLAST XML response
                align_xml_text = send_request(request,
                    pack_to_send, packs_at_all, fasta_hname+".fasta", tmp_fpath)

                # Get result tsv lines
                result_tsv_lines = parse_align_results_xml(align_xml_text,
                    packet["names"])

                seqs_processed += len( packet["names"] )

                # Write the result to tsv
                write_result(result_tsv_lines, tsv_res_path, acc_fpath, fasta_hname, outdir_path)
            #}
            pack_to_send += 1

            if seqs_processed >= probing_batch_size:#{
                remove_tmp_files(tmp_fpath)
                stop = True
                break
            #}
        #}
    #}
    remove_tmp_files(tmp_fpath)
    if stop:
        break
#}

print("\n {} sequences have been processed\n".format(seqs_processed))

print("Here are Genbank records that can be used for further sorting by 'barapost.py'.")
print("They are sorted by their occurence in probing batch:")

# Print accessions and record names sorted by occurence
# "-x[1][2]:": minus because we need descending order, [1] -- get tuple of "other information",
#   [2] -- get 3-rd element (occurence)
for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][2]):#{ 
    s_letter = "s" if other_info[2] > 1 else ""
    print(" {} sequence{} - {}, '{}'".format(other_info[2], s_letter, acc, other_info[1]))
#}

# Print number of unkmown sequences, if there are any:
unkn_num = probing_batch_size - sum( map(lambda x: x[2], acc_dict.values()) )
if unkn_num > 0:#{
    s_letter = "s" if unkn_num > 1 else ""
    print(" {} sequence{} - No significant similarity found".format(unkn_num, s_letter))
#}


print("""\nThey are saved in following file:
    '{}'""".format(acc_fpath))
print("""\nYou can edit this file before running 'barapost.py' in order to
  modify list of sequences that will be downloaded from Genbank
  and used as local (on your local computer) database used by 'barapost.py'.""")
end_time = time()
print('\n'+get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(end_time))) + "- Probing task is completed successfully!\n")
platf_depend_exit(0)