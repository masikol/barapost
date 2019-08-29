#!/usr/bin/env python3

# Version 2.1
# 29.08.2019 edition

# |===== Check python interpreter version =====|

from sys import version_info as verinf

if verinf.major < 3:#{
    print( "\nYour python interpreter version is " + "%d.%d" % (verinf.major, verinf.minor) )
    print("\tPlease, use Python 3!\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    raw_input("Press ENTER to exit:")
    exit(1)
#}

print("\n |=== barapost.py (version 2.1) ===|\n")

# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, sleep
start_time = time()


def get_work_time():#{
    return strftime("%H:%M:%S", localtime( time() ))
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
    "Function for printing error messages"
    print("\n   \a!! - ERROR: " + text + '\n')
#}


# |===== Handle command line arguments =====|
help_msg = """
DESCRIPTION:\n
 barapost.py -- script is designed for determinating the taxonomic position
    of nucleotide sequences by blasting each of them with 'blastn' from 'blast+' toolkit
    and regarding the best hit.\n
 'barapost.py' is meant to be used just after 'prober.py'.
 Script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.\n
 Results of the work of this script are written to TSV file,
    that can be found in result directory.\n
 FASTQ files processed by this script are meant to be sorted afterwards by 'fastQA_sorted.py'.
----------------------------------------------------------

Default parameters:\n
- all FASTQ and FASTA files in current directory will be processed;
- packet size (see '-p' option): 100 sequences;
- algorithm (see '-a' option): 'megaBlast';

  Default behavior is to download records-hits from Genbank according to results
of work of 'prober.py' script, build an indexed local database which consists of
downloaded sequences, and continue aligning with 'blast+' toolkit in order to save time.
----------------------------------------------------------

OPTIONS (* means mandatory option):\n
    -h (--help) --- show help message;\n
    -f (--infile) --- input FASTQ or FASTA (or '.fastq.gz', '.fasta.gz') file;
        You can specify multiple input files with this option (see EXAMPLES #2);\n
    -d (--indir) --- directory which contains FASTQ of FASTA files meant to be processed.
        E.i. all FASTQ and FASTA files in this direcory will be processed;
        Input files can be gzipped.\n
    -p (--packet-size) --- size of the packet, e.i. number of sequence to blast in one request.
        Value: integer number [1, 500]. Default value is 100;\n
    -a (--algorithm) --- BLASTn algorithm to use for aligning.
        Available values: 'megaBlast', 'discoMegablast', 'blastn'.
        Default is megaBlast;\n
  * -r (--prober-result-dir) --- directory with results of script 'prober.py'
        This is directory specified to 'prober.py' by '-o' option.
        Or named 'prober_result' if you've ran 'prober.py' with default options.
----------------------------------------------------------

EXAMPLES:\n
  1) Process one FASTQ file with default settings.
     File 'reads.fastq' has been already processed by 'prober.py'.
     Results of prober.py work are in directory 'prober_outdir':\n
       ./barapost.py -f reads.fastq -r prober_outdir\n
  2) Process FASTQ file and FASTA file with discoMegablast, packet size of 5 sequences.
     Files 'reads.fastq.gz' and 'another_sequences.fasta' have been already processed by 'prober.py'.
     Results of prober.py work are in directory 'prober_outdir':\n
       ./barapost.py -f reads.fastq.gz -f another_sequences.fasta -a discoMegablast -r prober_outdir\n
  3) Process all FASTQ and FASTA files in directory named "some_dir".
  All these files have been already processed by 'prober.py'.
  Results of prober.py work are in directory 'prober_outdir':\n
       ./barapost.py -d some_dir -r prober_outdir
"""
from sys import argv
import getopt

try:#{
    opts, args = getopt.getopt(argv[1:], "hf:d:o:p:a:r:",
        ["help", "infile=", "indir=", "packet-size=", "algorithm=", "prober-result-dir="])
#}
except getopt.GetoptError as gerr:#{
    print( str(gerr) )
    platf_depend_exit(2)
#}

is_fq_or_fa = lambda f: True if not re_search(r".*\.f(ast)?(a|q)(\.gz)?$", f) is None else False

# Default values:
fq_fa_list = list()
indir_path = None
packet_size = 100
blast_algorithm = "megaBlast"
prober_res_dir = None

if len(args) != 0:#{
    print_error("barapost.py does not take any positional arguments")
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

    if opt in ("-r", "--prober-result-dir"):#{
        if not os.path.exists(arg):#{
            print_error("directory '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        #}
        if not os.path.isdir(arg):#{
            print_error("'{}' is not a directory!".format(arg))
            platf_depend_exit(1)
        #}
        prober_res_dir = arg
    #}
#}

# Check if prober.py result directory is specified
if prober_res_dir is None:#{
    print("\n\t\a!! - ERROR: '-r' option is mandatory!\n")
    platf_depend_exit(1)
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


# Check if 'blast+' tookit is installed
pathdirs = os.environ["PATH"].split(os.pathsep)

for utility in ("blastn", "makeblastdb", "makembindex"):#{

    utility_found = False

    for directory in pathdirs:#{
        if os.path.exists(directory) and utility in os.listdir(directory):#{
            utility_found = True
            break
        #}
    #}

    if not utility_found:#{
        print("\tAttention!\n'{}' from blast+ toolkit is not found in your system.".format(utility))
        print("""If this error still occure although you have installed everything 
-- make sure that this program is added to PATH)""")
        platf_depend_exit(1)
    #}
#}


acc_fpath = os.path.join(prober_res_dir, "{}_probe_acc_list.tsv".format(blast_algorithm)) # form path to accession file

if not os.path.exists(acc_fpath):#{
    print_error("accession file '{}' not found!".format(acc_fpath))
    platf_depend_exit(1)
#}


print( strftime("\n%H:%M:%S", localtime(start_time)) + " - Start working\n")


# |===== Function for checking if 'https://ncbi.nlm.nih.gov' is available =====|

def check_connection(address):#{
    """
    Function checks if 'https://ncbi.nlm.nih.gov' is available.

    :return: None if 'https://ncbi.nlm.nih.gov' is available;
    """

    try:#{

        ncbi_server = address
        status_code = urllib.request.urlopen(ncbi_server).getcode()

        # Just in case
        if status_code != 200:#{
            print('\n' + get_work_time() + " - Site '{}' is not available.".format(address))
            print("Check your Internet connection.\a")
            print("Status code: {}".format(status_code))
            platf_depend_exit(-2)
        #}
        return
    #}
    except OSError as err:#{

        print('\n' + get_work_time() + " - Site '{}' is not available.".format(address))
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


# |===== Functionality for proper processing of gzipped files =====|

OPEN_FUNCS = (open, open_as_gzip)

is_gzipped = lambda file: True if file.endswith(".gz") else False

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


get_phred33 = lambda q_symb: ord(q_symb) - 33

def get_read_sum_qual(qual_str):#{

    phred33 = map(get_phred33, list(qual_str))
    sum_qual = sum(phred33)
    return sum_qual
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
    fasta_path = os.path.join(new_dpath, fasta_path) # place FASTA file into result directory

    fastq_patt = r".*\.f(ast)?q(\.gz)?$"

    num_lines = 0 # variable for counting lines in a file
    if not re_search(fastq_patt, fq_fa_path) is None and not os.path.exists(fasta_path):#{

        # Get ready to process gzipped files
        how_to_open = OPEN_FUNCS[ is_gzipped(fq_fa_path) ]
        fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_fa_path) ]

        with how_to_open(fq_fa_path) as fastq_file, open(fasta_path, 'w') as fasta_file:#{

            sum_qual = 0
            total_len = 0

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
                    sum_qual += get_read_sum_qual(line.strip())
                    total_len += len(line.strip())
                    counter = 0
                counter += 1
                num_lines += 1
            #}
        #}
        num_reads = int(num_lines / FASTQ_LINES_PER_READ) # get number of sequences
        file_avg_qual = round(sum_qual / total_len, 2)

        print("\n{}. '{}' ({} reads) --> FASTA".format(i+1, os.path.basename(fq_fa_path), num_reads))
        print("\tAverage quality of reads in this file: {} (Phred33)".format(file_avg_qual))
    #}
    # IF FASTA file is already created
    # We need only number of sequences in it.
    elif not re_search(fastq_patt, fq_fa_path) is None and os.path.exists(fasta_path):#{
        num_lines = sum(1 for line in open(fq_fa_path, 'r')) # get number of lines
        num_reads = int( num_lines / FASTQ_LINES_PER_READ ) # get number of sequences
    #}
    # If we've got FASTA source file
    # We need only number of sequences in it.
    else:#{
        num_lines = sum(1 for line in open(fq_fa_path, 'r')) # get number of lines
        num_reads = int( num_lines / FASTA_LINES_PER_SEQ ) # get number of sequences
        # If we've got FASTA source file we do not need to copy it
        fasta_path = fq_fa_path
    #}

    print("\n |========== file: '{}' ===========|".format(os.path.basename(fasta_path)))
    return {"fpath": fasta_path, "nreads": num_reads}
#}


def rename_file_verbosely(file, directory):#{
    is_analog = lambda f: file[file.rfind('.')] in f
    num_analog_files = len( list(filter(is_analog, os.listdir(directory))) )
    try:#{
        print('\n' + get_work_time() + " - Renaming old file:")
        name_itself = file[: file.rfind('.')]
        ext = file[file.rfind('.'):]
        num_analog_files = str(num_analog_files)
        new_name = name_itself+"_old_"+num_analog_files+ext
        print("\t'{}' --> '{}'".format(os.path.basename(file), new_name))
        os.rename(file, new_name)
    #}
    except:#{
        # Anything (and not only strings) can be passed to the function
        print("\nFile '{}' cannot be renamed".format( str(file)) )
    #}
#}


def look_around(new_dpath, fasta_path, blast_algorithm):#{
    """
    Function looks around in order to ckeck if there are results from previous runs of this script.

    Returns None if there is no result from previous run.
    If there are results from previous run, returns a dict of the following structure:
    {
        "pack_size": packet_size (int),
        "attmpt": saved_attempt (int),
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

    # "hname" means human readable name (e.i. without file path and extention)
    fasta_hname = os.path.basename(fasta_path) # get rid of absolute path
    fasta_hname = fasta_hname[: fasta_hname.rfind(".fasta")] # get rid of '.fasta' extention

    # Form path to temporary file
    tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)
    # Form path to result file
    tsv_res_fpath = "{}_{}_result.tsv".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)

    print("\nWriting results to file '{}'\n".format(tsv_res_fpath))

    num_done_reads = None # variable to keep number of succeffdully processed sequences

    if os.path.exists(tsv_res_fpath):#{
        with open(tsv_res_fpath, 'r') as res_file:#{
            try:#{ There can be invalid information in result file
                lines = res_file.readlines()
                num_done_reads = len(lines) - 1 # the first line is a head
                last_line = lines[-1]
                last_seq_id = last_line.split(DELIM)[0]

                # # Get accessions from previous run
                # global acc_dict
                # for line in lines[1:]:#{ omit table head
                #     acc = line.strip().split(DELIM)[2]
                #     acc_dict[acc] = line.strip().split(DELIM)[1]
                # #}
            #}
            except Exception as err:#{
                print("\nData in result file '{}' is broken.".format(tsv_res_fpath))
                print( str(err) )
                print("Start from the beginning.")
                rename_file_verbosely(tsv_res_fpath, new_dpath)
                rename_file_verbosely(tmp_fpath, new_dpath)
                return None
            #}
            else:#{
                print("Last processed sequence: " + last_seq_id)
            #}
        #}
    #}

    # If we start from the beginning, we have no sequences processed
    if num_done_reads is None:
        num_done_reads = 0

    if os.path.exists(tmp_fpath):#{
        try:#{ There can be invalid information in tmp file of tmp file may not exist
            with open(tmp_fpath, 'r') as tmp_file:
                temp_lines = tmp_file.readlines()
            packet_size = int(temp_lines[0])
            print("\nPacket size switched to '{}', as it was during previous run".format(packet_size))
            attempt_save = int(temp_lines[1].split(DELIM)[0])
        #}
        except Exception as err:#{
            print("\nData in temporary file '{}' is broken.".format(tmp_fpath))
            print( str(err) )
            print("Start from the beginning.")
            rename_file_verbosely(tsv_res_fpath, new_dpath)
            rename_file_verbosely(tmp_fpath, new_dpath)
            return None
        else:#{
            # Return data from previous run
            return {
                "pack_size": packet_size,
                "attmpt": attempt_save,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        #}
    #}
    else:#{
        return
    #}
#}


def get_packet(fasta_file, packet_size):#{
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

    packet = ""
    names = list()

    for i in range(packet_size):#{

        seq_id = fasta_file.readline().strip().partition(' ')[0] # prune the seq id a little bit
        seq = fasta_file.readline().strip() # get the sequence

        # if 'seq_d is an empty string -- 'seq' is as well
        if seq_id == "":  # if previous sequence was last in current file
            break
        else:
            names.append(seq_id)
            packet += "{}\n{}\n".format(seq_id, seq)
    #}

    return {"fasta": packet.strip(), "names": names}  # remove the last '\n' character
#}

# Variable for counting accessions of records menat to be downloaded from Genbank.
# Is used only for printing the list of accessions to console.
acc_counter = 0

import threading # There is no need to import this module if '--remote-only' is specified.

def get_gi_by_acc(acc):#{
    """
    Function accepts sequence accession and returns it's GI number for further downloading
        with E-utilities.

    :param acc: sequence accession;
    :type acc: str;

    Returns sequence's GI number of 'str'.
    """
    
    server = "www.ncbi.nlm.nih.gov"
    url = "/nuccore/{}?report=gilist&log$=seqview&format=text".format(acc)

    error = True
    while error:#{
        try:#{
            conn = http.client.HTTPSConnection(server)
            conn.request("GET", url)
            response = conn.getresponse()
            gi_text = str(response.read(), "utf-8")
            # The only text in this HTML document is GI. It is in <pre> tag.
            gi = str( re_search(r"<pre>([0-9]+)", gi_text).group(1) )
            conn.close()
        #}
        except OSError as oserr:#{
            # If site is unavailable -- try to connect again in 60 seconds
            print_error("unable to connect to 'www.ncbi.nlm.nih.gov'")
            print("  barapost will try to connect adaing in 60 seconds")
            sleep(60)
            continue
        #}
        else:#{
            error = False
        #}
    #}

    # Print accession and record name to console
    global acc_counter
    acc_counter += 1
    print("\r {}/{} GI numbers are got".format( acc_counter, len(acc_dict.keys()) ), end="")

    return gi
#}


def retrieve_fastas_by_gi(gi_list, db_dir):#{
    """
    Function downloads set of records from Genbank according to list of GIs passed to it.
    Downloaded FASTA file will be placed in 'db_dir' directory and named 'local_seq_set.fasta'

    :param gi_list: list of GI numbers of sequences meant to be downloaded;
    :type gi_list: list<str>;
    :param db_dir: path to directory in which downloaded FASTA file will be placed;
    :type db_dir: str;

    Returns path to downloaded FASTA file of 'str'.
    """

    local_fasta = os.path.join(db_dir, "local_seq_set.fasta") # path to downloaded FASTA file
    gis_del_comma = ','.join(gi_list) # GI numbers must be separated by comma in url
    # E-utilities provide us with possibility of downloading records from Genbank by GI numbers.
    retrieve_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text".format(gis_del_comma)
    
    global stop_wait # a flag variable that will signal waiter-function to stop executing

    def download_waiter():#{
        """
        Function waits untill 'local_fasta' file is downloaded.
        It prints size of downloaded data to console during downloading.
        This function just waits -- it won't bring you the menu :).
        """
        # Wait untill downloading starts
        while not os.path.exists(local_fasta):
            if stop_wait:
                return
            sleep(1)
        print()

        fsize = round(os.path.getsize(local_fasta) / (1024**2), 1) # get megabytes
        while not stop_wait:#{
            # Get size of downloaded data
            fsize = round(os.path.getsize(local_fasta) / (1024**2), 1)
            print("\r {} MB downloaded".format(fsize), end="")
            sleep(1) # instant updates are not necessary
        #}
        # Print total size of downloaded file
        fsize = round(os.path.getsize(local_fasta) / (1024**2), 1)
        print("\r {} MB downloaded \n".format(fsize))
    #}

    waiter = threading.Thread(target=download_waiter) # create thread
    stop_wait = False # raise the flag
    waiter.start() # start waiting
    
    print("Downloading sequences for local database building...")

    urllib.request.urlretrieve(retrieve_url, local_fasta) # retrieve FASTA file
    stop_wait = True # lower the flag
    waiter.join() # main thread will wait until waiter function ends it's work

    print("Downloading is completed\n\n")

    return local_fasta
#}


def build_local_db(acc_dict, prober_res_dir):#{
    """
    Function builds a local indexed database with utilities from 'blast+' toolkit.

    :param acc_dict: a dictionary of accessions and record names
        Accession are keys, record names are values;
    :type acc_dict: dict<str, str>;
    :param prober_res_dir: path to current result directory (each processed file has it's own result directory);
    :type prober_res_dir: str;

    Returns path to builded local indexed database.
    """

    print("""\nFollowing sequences will be downloaded from Genbank
for further blasting on your local machine with 'blast+' toolkit:\n""")
    for i, acc in enumerate(acc_dict.keys()):#{
        print(" {}. {} - '{}'".format(i+1, acc, acc_dict[acc]))
    #}
    print() # just print blank line
    
    print('~'*20 + '\n')

    db_dir = os.path.join(prober_res_dir, "local_database") # path to directory in which database will be placed
    try:#{
        os.makedirs(db_dir)
    #}
    except OSError as err:#{
        #If this directory exists

        while True:#{

            print("Database directory exists.")
            if len(os.listdir(db_dir)) == 0:#{
                # If it is empty -- nothing stops us. break and build a database
                print("It is empty. Building a database...")
                break
            #}
            else:#{
                # If there are some files, the user will decide, whether to use this database
                # or to remove it and to build again (e.g. if the database haven't been builded successfully) .
                print("Here are files located in this directory:")
                for i, file in enumerate(os.listdir(db_dir)):
                    print("  {}. '{}'".format(i+1, file))

                reply = input("""\nPress ENTER to continue aligning using this database.
            Enter 'r' to rebuild database:>>""")

                if reply == "":
                    # Do not build a database, just return path to it.
                    return os.path.join(db_dir, "local_seq_set.fasta")
                elif reply == 'r':
                    # Empty this directory and break from the loop in order to build a database.
                    for file in os.listdir(db_dir):
                        os.unlink(file)
                    break
                else:
                    # Ask again
                    continue
            #}
        #}

    #}

    # Get list of GI numbers. Function 'get_gi_by_acc' will print the list of GIs to console.
    print("Getting GI numbers necessary for downloading sequences from Genbank...")
    gi_list = list( map(get_gi_by_acc, acc_dict.keys()) )
    print()

    local_fasta = retrieve_fastas_by_gi(gi_list, db_dir) # download FASTA file

    print("Creating database...")
    # Configure command line
    make_db_cmd = "makeblastdb -in {} -parse_seqids -dbtype nucl".format(local_fasta)
    exit_code = os.system(make_db_cmd) # make a blast-format database
    if exit_code != 0:#{
        print_error("error while making the database")
        platf_depend_exit(exit_code)
    #}
    print("""Database is successfully installed:
'{}'\n""".format(local_fasta))

    print("Creating database index...")
    # Configure command line
    make_index_cmd = "makembindex -input {} -iformat blastdb -verbosity verbose".format(local_fasta)
    exit_code = os.system(make_index_cmd) # create an index for the database
    if exit_code != 0:#{
        print_error("error while creating database index")
        platf_depend_exit(exit_code)
    #}
    print("Database index has been successfully created\n")

    # Gzip downloaded FASTA file in order to save space on disk
    with open(local_fasta, 'r') as fasta_file, open_as_gzip(local_fasta+".gz", "wb") as fagz_file:#{
        fagz_file.write(bytes(fasta_file.readline(), "utf-8"))
    #}
    os.unlink(local_fasta) # remove source FASTA file, not the database

    return local_fasta
#}


def configure_request(packet, blast_algorithm):#{
    """
    Function meant to replace 'configure_request' one that works with BLAST server.
    """

    # Algorithms in 'blast+' are named in a little different way comparing to BLAST server.
    if blast_algorithm == "megaBlast":
        blast_algorithm = "megablast"
    elif blast_algorithm == "discoMegablast":
        blast_algorithm = "dc-megablast"

    # Indexed discontiguous searches are not supported:
    # https://www.ncbi.nlm.nih.gov/books/NBK279668/#usermanual.Megablast_indexed_searches
    if blast_algorithm != "dc-megablast":
    	use_index = "true"
    else:
    	use_index = "false"

    # Confirue command line
    blast_cmd = "blastn -query {} -db {} -outfmt 5 -out {} -task {} -max_target_seqs 5 -use_index {}".format(query_path,
        local_fasta, xml_output_path, blast_algorithm, use_index)

    # Write packet sequences to query file
    with open(query_path, 'w') as query_file:#{
        query_file.write(packet)
    #}

    return blast_cmd # return command line in order to follow the interface used in 'kernel loop'
#}

def send_request(request, attempt, attempt_all, filename, tmp_fpath):#{
    """
    Function meant to replace 'send_request' one that works with BLAST server.
    """

    # Save temporary data
    with open(tmp_fpath, 'a') as tmpfile:
        tmpfile.write("{}\n".format(attempt))

    print("\n{} - Alignment for: '{}' ({}/{}) started".format(get_work_time(), filename, attempt, attempt_all))

    # request is blastn cmd returned by 'configure_localdb_request' function
    exit_code = os.system(request)
    if exit_code != 0:#{
        print_error("error while aligning a sequence against local database")
        platf_depend_exit(exit_code)
    #}

    # Read response XML text from file
    with open(xml_output_path, 'r') as xml_file:#{
        align_xml_text = xml_file.read().strip()
    #}

    return align_xml_text # return it
#}


query_path = "query.fasta"  # path to file in which 'blastn will write it's report
xml_output_path = "xml_output.xml"  # path to file in which a packet of sequences will be written

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

    # Iterate through "Iteration" and "Iteration_hits" nodes
    for iter_elem, iter_hit in zip(root.iter("Iteration"), root.iter("Iteration_hits")):#{
    
        # "Iteration" node contains query name information
        query_name = iter_elem.find("Iteration_query-def").text

        query_len = iter_elem.find("Iteration_query-len").text

        # If there are any hits, node "Iteration_hits" contains at least one "Hit" child
        hit = iter_hit.find("Hit")
        if hit is not None:#{

            # Get full hit name (e.g. "Erwinia amylovora strain S59/5, complete genome")
            hit_name = hit.find("Hit_def").text
            # Format hit name (get rid of stuff after comma)
            hit_taxa_name = hit_name[: hit_name.find(',')] if ',' in hit_name else hit_name
            hit_taxa_name = hit_taxa_name.replace(" complete genome", "") # sometimes there are no comma before it
            hit_taxa_name = hit_taxa_name.replace(' ', '_')


            hit_acc = hit.find("Hit_accession").text # get hit accession

            # Find the first HSP (we need only the first one)
            hsp = next(hit.find("Hit_hsps").iter("Hsp"))

            align_len = hsp.find("Hsp_align-len").text.strip()

            pident = hsp.find("Hsp_identity").text # get number of matched nucleotides

            gaps = hsp.find("Hsp_gaps").text # get number of gaps

            evalue = hsp.find("Hsp_evalue").text # get e-value
            # If E-value is low enough -- add this subject sequence to 'acc_dict' to further downloading

            pident_ratio = round( float(pident) / int(align_len) * 100, 2)
            gaps_ratio = round( float(gaps) / int(align_len) * 100, 2)

            print("""\n '{}' -- '{}'
    Query length - {} nt;
    Identity - {}/{} ({}%); Gaps - {}/{} ({}%)""".format(query_name, hit_taxa_name,
                    query_len, pident, align_len, pident_ratio, gaps, align_len, gaps_ratio))

            # Append new tsv line containing recently collected information
            result_tsv_lines.append( DELIM.join( [query_name, hit_taxa_name, hit_acc, query_len,
                align_len, pident, gaps, evalue] ))
        #}
        else:
            # If there is no hit for current sequence
            print("\n {} -- No significant similarity found.\n Query length - {}.".format(query_name, query_len))
            result_tsv_lines.append(DELIM.join( [query_name, "No significant similarity found", "", query_len, ] ))
    #}
    return result_tsv_lines
#}


def write_result(res_tsv_lines, tsv_res_path):#{
    """
    Function writes result of blasting to result tsv file.

    :param res_tsv_lines: tsv lines returned by 'parse_align_results_xml()' funciton;
    :type res_tsv_lines: list<str>;
    :param tsv_res_path: path to reslut tsv file;
    :type tsv_res_path: str;
    """


    # If there is no result tsv fil -- create it and write a head of the table.
    if not os.path.exists(tsv_res_path):#{
        with open(tsv_res_path, 'w') as tsv_res_file:#{
            tsv_res_file.write(DELIM.join( ["QUERY_ID", "HIT_NAME", "HIT_ACCESSION", "QUERY_LENGTH",
                "ALIGNMENET_LENGTH", "IDENTITY", "GAPS", "E-VALUE"] ) + '\n')
        #}
    #}
    # Write reslut tsv lines to this file
    with open(tsv_res_path, 'a') as tsv_res_file:#{
        for line in res_tsv_lines:
            tsv_res_file.write(line + '\n')
    #}
#}

def get_curr_res_dir(fq_fa_path, prober_res_dir):#{
    
    # dpath means "directory path"
    new_dpath = os.path.join(prober_res_dir, os.path.basename(fq_fa_path)) # get rid of absolute path
    new_dpath = new_dpath[: new_dpath.rfind(".fast")] # get rid of extention

    return new_dpath
#}

def configure_acc_dict(acc_fpath):#{

    acc_dict = dict()

    with open(acc_fpath, 'r') as acc_file:#{
        lines = acc_file.readlines()

        for line in lines:#{
            if line.strip() != "" and not line.strip().startswith('#') and not line.startswith("ACCESSION"):#{
                acc = intern(line.strip().split(DELIM)[0])
                name = line.strip().split(DELIM)[1]
                acc_dict[acc] = name
            #}
        #}
    #}

    if len(acc_dict) == 0:#{
        print_error("no accession information found in file '{}".format(acc_fpath))
        platf_depend_exit(1)
    #}

    return acc_dict
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
#        "attmpt": saved_attempt (int),
#        "tsv_respath": path_to_tsv_file_from_previous_run (str),
#        "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file (int),
#        "tmp_fpath": path_to_pemporary_file (str)
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. 'packet' is a dict of the following structure:
#    {
#        "fasta": FASTA_data_containing_query_sequences (str),
#        "names": list_of_sequence_ids_from_FASTA_file (list<str>)
#    }

#                   |===== Kernel loop =====|

print("\nChecking Internet connection...")
check_connection("https://ncbi.nlm.nih.gov")
print("OK\n")

print(" Following files will be processed:")
for i, path in enumerate(fq_fa_list):#{
    print("    {}. '{}'".format(i+1, path))
#}
print('-'*30 + '\n')

# It is a dictionary of accessions and record names.
# Accessions are keys, record names are values.
# This dictionary is filled while processing and at the beginning of continuation.
acc_dict = configure_acc_dict(acc_fpath)

# Build a database
local_fasta = build_local_db(acc_dict, prober_res_dir)

# Iterate through found source FASTQ and FASTA files
for i, fq_fa_path in enumerate(fq_fa_list):#{
    
    # Create the result directory with the name of FASTQ of FASTA file being processed:
    new_dpath = get_curr_res_dir(fq_fa_path, prober_res_dir)

    # Convert FASTQ file to FASTA (if it is FASTQ) and get it's path and number of sequences in it:
    curr_fasta = fastq2fasta(fq_fa_path, i, new_dpath)

    # "hname" means human readable name (e.i. without file path and extention)
    fasta_hname = os.path.basename(curr_fasta["fpath"]) # get rid of absolure path
    fasta_hname = fasta_hname[: fasta_hname.rfind(".fasta")] # get rid of file extention

    # Look around and ckeck if there are results of previous runs of this script
    # If 'look_around' is None -- there is no data from previous run
    previous_data = look_around(new_dpath, curr_fasta["fpath"],
        blast_algorithm)

    if previous_data is None:#{ # If there is no data from previous run

        num_done_reads = 0 # number of successfully processed sequences
        saved_attempt = None # number of last successful attempt (there is no such stuff for de novo run)
        tsv_res_path = "{}_{}_result.tsv".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form result tsv file path
        tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form temporary file path
    #}
    else:#{ # if there is data from previous run

        num_done_reads = previous_data["n_done_reads"] # get number of successfully processed sequences
        saved_attempt = previous_data["attmpt"] # get number of last successful attempt
        packet_size = previous_data["pack_size"] # packet size sholud be the same as it was in previous run
        tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        tmp_fpath = previous_data["tmp_fpath"] # temporary file sholud be the same as during previous run
    #}

    attempt_all = curr_fasta["nreads"] // packet_size # Calculate total number of packets sent from current FASTA file
    if curr_fasta["nreads"] % packet_size > 0: # And this is ceiling (in order not to import 'math')
        attempt_all += 1
    attempts_done = int( num_done_reads / packet_size ) # number of successfully processed sequences

    with open(curr_fasta["fpath"], 'r') as fasta_file:#{

        # Go untill the last processed sequence
        for _ in range( int(num_done_reads * 2) ):#{
            fasta_file.readline()
        #}

        reads_left = curr_fasta["nreads"] - num_done_reads # number of sequences left to precess
        attempts_left = attempt_all - attempts_done # number of packets left to send
        attempt = attempts_done+1 if attempts_done > 0 else 1 # current attempt

        # Iterate througth packets left to send
        for i in range(attempts_left):#{

            with open(tmp_fpath, 'w') as tmp_file:
                tmp_file.write(str(packet_size)+ '\n')

            packet = get_packet(fasta_file, packet_size) # form the packet

            if packet["fasta"] is "":#{   Just in case
                print("Recent packet is empty")
                break
            #}

            print("\nGo to BLAST (" + blast_algorithm + ")!")
            print("Request number {} out of {}.".format(attempt, attempt_all))
            
            request = configure_request(packet["fasta"], blast_algorithm) # get the request

            # Send the request get BLAST XML response
            align_xml_text = send_request(request,
                attempt, attempt_all, fasta_hname+".fasta", tmp_fpath)

            # Get result tsv lines
            result_tsv_lines = parse_align_results_xml(align_xml_text,
                packet["names"])

            # Write the result to tsv
            write_result(result_tsv_lines, tsv_res_path)

            attempt += 1
        #}
    #}
    remove_tmp_files(tmp_fpath)
#}
remove_tmp_files(query_path, xml_output_path)


print("\nTask completed successfully!")
platf_depend_exit(0)