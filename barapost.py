#!/usr/bin/env python3

# Version 1.0
# 22.08.2019 edition

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

print("\n |=== barapost.py (version 1.0) ===|\n")

# |===== Function for dealing with time =====|

from time import time, strftime, localtime, sleep
start_time = time()

print( strftime("\n%H:%M:%S", localtime(start_time)) + " - START WORKING\n")

def get_work_time():#{
    return strftime("%H:%M:%S", localtime( time() ))
#}

# |===========================================|

import os
from re import search as re_search
from gzip import open as open_as_gzip # input files might be gzipped
from xml.etree import ElementTree

import http.client
import urllib.request
from urllib.error import HTTPError
import urllib.parse


# |===== Function that asks to press ENTER on Windows =====|

from sys import platform

def platf_depend_exit(exit_code):#{
    """
    A function that asks to press ENTER press on Windows
        and exits.

    :type exit_code: int;
    """
    if platform.startswith("win"):
        input("Press ENTER to exit:")
    exit(exit_code)
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
        print( repr(err) )

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
                print("\n\tNot a VALID number entered!\a\n" + '~'*20)
                continuation = None
            #}
            else:#{
                print("You have chosen number " + str(continuation) + '\n')
                print('~' * 20 + '\n')
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
    limit = num_reads if num_reads <= 1000 else 1000

    while packet_size is None:#{
        
        packet_size = input("""
Please, specify the number of sequences that should be sent to the NCBI server in one request.
E.g. if you have 10 reads in your file, you can send 10 reads as single
    request -- in this case you should enter number 10. You may send 2 requests containing
    5 reads both -- in this case you should enter number 5.


There are {} reads left to process in current file.
Enter the number (from 1 to {}):>> """.format(num_reads, limit))
        # Check if entered value is integer number. If no, give another attempt.
        try:#{
            packet_size = int(packet_size)
            if packet_size < 1 or packet_size > limit:#{ Check if input number is in [1, limit]
                print("\n\tNot a VALID number entered!\a\n" + '~'*20)
                packet_size = None
            #}
            else:#{
                print("You have chosen number " + str(packet_size) + '\n')
                print('~' * 20 + '\n')
            #}
        #}
        except ValueError:#{
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            packet_size = None
        #}
    #}
    return(packet_size)
#}


def get_algorithm():#{
    """
    Function asks the user about what BLASTn algorithm to use.

    :return: one of the following strings "megaBlast", "discoMegablast", "blastn";
    :return type: str;
    """

    reply = None
    blast_algorithm = "megaBlast" # default value

    while reply is None:#{
        reply = input("""
Please choose a BLAST algorithm:
    1. Highly similar sequences (megablast)
    2. Optimize for More dissimilar sequences (discontiguous megablast)
    3. Optimize for Somewhat similar sequences (blastn)

Enter the number (1, 2 or 3):>> """)
        # Check if entered value is integer number. If no, give another attempt.
        try:#{
            reply = int(reply)
            if reply < 1 or reply > 3:#{ Check if input number is in [1, 3]
                print("\n\tNot a VALID number entered!\a\n" + '~'*20)
                reply = None
            #}
            else:#{
                print("You have chosen number "+ str(reply) + '\n')
                print('~' * 20 + '\n')

                if reply == 1:
                    blast_algorithm = "megaBlast"
                elif reply == 2:
                    blast_algorithm = "discoMegablast"
                else:
                    blast_algorithm = "blastn"
            #}
        #}
        except ValueError:
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            reply = None
            algorithm = None
    #}
    return blast_algorithm
#}

def get_organisms():#{
    
    # Only 2 'nt' database restrictions for now
    max_org = 2
    organisms = list()
    
    error = True
    while error:#{
        reply = input("""\n~~~~~~~~~~~~~~~~~~~~~~~~~~~
Please specify organisms for 'nt' database restriction.
You are to enter organism's name and it's Taxonomy ID (aka taxid) separated by comma
 after prompt symbol '>>' as follows:\n
    >> Erwinia,551\n
You can enter 2 organisms separated by semicolon as follows:\n
    >> Erwinia,551;Pseudomonas,286

You can find taxid of your organism from NCBI taxonomy browser:
\t'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi'

Default settings are:    bacteria,2;viruses,10239

Enter organism's name and taxid:>> """)

        if reply is "":#{
            organisms.append("bacteria (taxid:2)") # default
            organisms.append("viruses (taxid:10239)") # default
            print("\nYou have chosen following organisms:")
            for i, org in enumerate(organisms):#{
                print("\t{}. {}".format(i+1, org))
            #}
            print('~'*30 + '\n')
            return organisms
        #}

        try:#{
            org_list = reply.strip().split(';')
            org_list = list( map(str.strip, org_list) )

            if len(org_list) > max_org:#{
                raise Exception("\nYou can enter only 1 or 2 organisms and they should be separated by semicolon (;)\a")
            #}

            for org in org_list:#{
                name_and_taxid = org.strip().split(',')
                name_and_taxid = list( map(str.strip, name_and_taxid) )
                if len(name_and_taxid) != 2:
                    raise Exception("\nOrganism's name and it's taxid should be separated by comma (,)\a")
                # Validate TaxID inteder format: it will raise ValueError if taxid is invalid
                tmp_taxid = int(name_and_taxid[1])
                if tmp_taxid < 1:
                    raise Exception("\nTaxID should be positive integer number\a")
                organisms.append("{} (taxid:{})".format(name_and_taxid[0], name_and_taxid[1]))
            #}
        #}
        except ValueError:#{
            print("\n" + "=/"*20 + "\n")
            print("\n\t!! - ERROR - TaxID should be positive integer number\a")
            sleep(2)
            continue
        #}
        except Exception as err:#{
            print("\n" + "=/"*20 + "\n")
            print("\n\t!! - ERROR: invalid input format")
            print( str(err) )
            sleep(2)
            continue
        #}

        print("\nYou have chosen following organisms:")
        for i, tax_fmt in enumerate(organisms):#{
            print("\t{}. {}".format(i+1, tax_fmt))
        #}
        print('~'*20 + '\n')
        error = False

    return organisms
#}


# |===== End of question funtions =====|


is_fastq_or_fastqgz = lambda f: f.endswith(".fastq") | f.endswith(".fastq.gz")


def get_fastq_list():#{

    # Get all '.fastq' and '.fastq.gz' files in current directory
    fastq_list = os.listdir('.')
    global is_fastq_or_fastqgz
    fastq_list = list( filter(is_fastq_or_fastqgz, fastq_list) )

    # If there are files to process
    if len(fastq_list) != 0:#{
        print("\nFollowing files have been found and will be processed:")
        for i, line in enumerate(fastq_list):#{
            print("\t{}. '{}'".format( str(i+1), fastq_list[i] ))
        #}
        print('-' * 30)
    #}
    else:#{
        print("\nNo '.fastq' or '.fastq.gz' files have been found in current directory!")
        platf_depend_exit(1) # exit
    #}
    return fastq_list
#}


# |===== Functionality for proper processing of gzipped files =====|

OPEN_FUNCS = (open, open_as_gzip)

is_gzipped = lambda file: True if file.endswith(".gz") else False

# Data from .fastq and .fastq.gz should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format '.fastq' line
    lambda line: line.decode("utf-8").strip()  # format '.fastq.gz' line
)

# |=== Delimiter for result tsv file ===|
DELIM = '\t'

# |=== File format constants ===|
FASTQ_LINES_PER_READ = 4
FASTA_LINES_PER_READ = 2


def fastq2fasta(fastq_path, i, new_dpath):#{
    """
    Function conwerts FASTQ file to FASTA format, if there is no FASTA file with
        the same name as FASTQ file.

    :param fastq_path: path to FASTQ file being processed;
    :type fastq_path: str;
    :param i: order number of fastq_file;
    :type i: int;
    :param new_dpath: path to current (corresponding to fastq_path file) result directory;
    :type new_dpath: str;

    Returns dict of the following structure:
    {
        "fpath": path_to_FASTA_file (str),
        "nreads": number_of_reads_in_this_FASTA_file (int)
    }
    """
    
    fasta_path = fastq_path.replace(".fastq", ".fasta") # change extention
    fasta_path = os.path.join(new_dpath, fasta_path) # place FASTA file into result directory

    global FASTQ_LINES_PER_READ
    global FASTA_LINES_PER_READ

    num_lines = 0 # variable for counting lines in a file
    if not os.path.exists(fasta_path):#{

        # Get ready to process gzipped files
        how_to_open = OPEN_FUNCS[ is_gzipped(fastq_path) ]
        fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq_path) ]

        with how_to_open(fastq_path) as fastq_file, open(fasta_path, 'w') as fasta_file:#{

            counter = 1 # variable for retrieving only 1-st and 2-nd line of FASTQ record
            for line in fastq_file:#{
                line = fmt_func(line)
                if counter <= 2:#{      write only 1-st and 2-nd line out of 4
                    if line[0] == '@':
                        line = '>' + line[1:]  # replace '>' with '@'
                    fasta_file.write(line + '\n')
                #}
                # reset the counter if the 4-th (quality) line has been encountered
                elif counter == 4:
                    counter = 0
                counter += 1
                num_lines += 1
            #}
        #}
        num_reads = int(num_lines / FASTQ_LINES_PER_READ) # get number of reads
    #}
    # If there if corresponding FASTA file from previous run, we do not need to create it.
    # We need only number of reads in it.
    else:#{
        num_lines = sum(1 for line in open(fasta_path, 'r')) # get number of lines
        num_reads = int( num_lines / FASTA_LINES_PER_READ ) # get number of reads
    #}

    print("\n{}. '{}' ({} reads) --> FASTA\n".format(i+1, fastq_path, num_reads))

    return {"fpath": fasta_path, "nreads": num_reads}
#}


def rename_file_verbosely(file, directory):#{
    is_analog = lambda f: file[file.rfind('.')] in f
    num_analog_files = len( list(filter(is_analog, os.listdir(directory))) )
    print('~' * 20)
    try:#{
        print('\n' + get_work_time() + " - Renaming old file:")
        name_itself = os.path.basename(file)[: os.path.basename(file).rfind('.')]
        ext = file[file.rfind('.'):]
        num_analog_files = str(num_analog_files)
        new_name = name_itself+"_old_"+num_analog_files+ext
        print("\t'{}' --> '{}'".format(file, new_name))
        os.rename(file, new_name)
    #}
    except:#{
        # Anything (and not only strings) can be passed to the function
        print("\nFile '{}' cannot be renamed".format( str(file)) )
    #}
    print('~' * 20 + '\n')
#}


def look_around(new_dpath, fasta_path, blast_algorithm):#{
    """
    Function looks around in order to ckeck if there are results from previous runs of this script.

    Returns None if there is no result from previous run.
    If there are results from previous run, returns a dict of the following structure:
    {
        "pack_size": packet_size (int),
        "attmpt": saved_attempt (int),
        "RID": saved_RID (str),
        "tsv_respath": path_to_tsv_file_from_previous_run (str),
        "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file (int),
        "tmp_fpath": path_to_pemporary_file (str)
    }

    :param new_dpath: path to current (corresponding to fastq_path file) result directory;
    :type new_dpath: str;
    :param fasta_path: path to current (corresponding to fastq_path file) FASTA file;
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

    num_done_reads = None # variable to keep number of succeffdully processed reads

    print("\n |========== file: '{}' ===========|".format(fasta_path))
    continuation = False    # default value
    # Check if there are results from previous run.
    if os.path.exists(tsv_res_fpath):#{
        print('\n' + get_work_time() + " - The previous result file is found in the directory:")
        print("\t'{}'".format(tmp_fpath))
        continuation = is_continued() # Allow politely to continue from last successful attempt.
        if not continuation:
            rename_file_verbosely(tsv_res_fpath, new_dpath)
            rename_file_verbosely(tmp_fpath, new_dpath)
    #}

    if continuation:#{   Find the name of last successfulli processed sequence
        print("Let's try to continue...")
        if os.path.exists(tsv_res_fpath):#{
            with open(tsv_res_fpath, 'r') as res_file:#{
                try:#{ There can be invalid information in result file
                    lines = res_file.readlines()
                    num_done_reads = len(lines) - 1 # the first line is a head
                    last_line = lines[-1]
                    last_seq_id = last_line.split(DELIM)[0]
                #}
                except Exception as err:#{
                    print("\nData in result file '{}' is broken.".format(tsv_res_fpath))
                    print("Start from the beginning.")
                    rename_file_verbosely(tsv_res_fpath, new_dpath)
                    rename_file_verbosely(tmp_fpath, new_dpath)
                    return None
                #}
                else:#{
                    print("Last successful attempt: " + last_seq_id)
                #}
            #}
        #}

        # If we start from the beginning, we have no reads processed
        if num_done_reads is None:
            num_done_reads = 0

        # Get packet size, number of the last attempt and RID
        try:#{ There can be invalid information in tmp file of tmp file may not exist
            with open(tmp_fpath, 'r') as tmp_file:
                temp_lines = tmp_file.readlines()
            packet_size = int(temp_lines[0])
            attempt_save = int(temp_lines[-1].split(DELIM)[0])
            RID_save = temp_lines[-1].split(DELIM)[1].strip()
        #}
        except Exception:#{
            print("\nTemporary file '{}' not found of broken!".format(tmp_fpath))
            print("{} reads have been already processed".format(num_done_reads))
            total_num_seqs = int( sum(1 for line in open(fasta_path, 'r')) / FASTA_LINES_PER_READ )
            print("{} reads left".format(total_num_seqs - num_done_reads))
            packet_size = get_packet_size(total_num_seqs - num_done_reads)
            return {
                "pack_size": packet_size,
                "attmpt": int(num_done_reads / packet_size),
                "RID": None,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        else:#{
            # Return data from previous run
            return {
                "pack_size": packet_size,
                "attmpt": attempt_save,
                "RID": RID_save,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_reads,
                "tmp_fpath": tmp_fpath
            }
        #}
        #}

    # Remove files from previous run if we've desided to start from the beginning
    else:#{
        return None
    #}
#}


def get_packet(fasta_file, packet_size):#{
    """
    Function collects the packet of query sequences to send to BLAST server.

    :param fasta_file: file instance of FASTA file to retrieve reads from;
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


def configure_request(packet, blast_algorithm, organisms):#{
    """
    Function configures the request to BLAST server.

    :param packet: FASTA_data_containing_query_sequences;
    :type packet: str;
    :param blast_algorithm: BLASTn algorithm to use;
    :type blast_algorithm: str;
    :param organisms: list of strings performing 'nt' database restriction;
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

    # 'nt' database restrictions:
    for i, org in enumerate(organisms):#{
        payload["EQ_MENU{}".format(i if i > 0 else "")] = org
    #}
    payload["NUM_ORG"] = str( len(organisms) )

    payload = urllib.parse.urlencode(payload)
    headers = { "Content-Type" : "application/x-www-form-urlencoded" }

    return {"payload":payload, "headers": headers}
#}

def send_request(request):#{
    """
    Function sends a request to "blast.ncbi.nlm.nih.gov/blast/Blast.cgi".

    :param request: request_data (it is a dict that 'configure_request()' function returns);
    :param request: dict<dict>;

    Returns a dict of the following structure:
    {
        "RID": Request ID (str),
        "RTOE", time_to_wait_provided_by_the_NCBI_server (int)
    }
    """
    payload = request["payload"]
    headers = request["headers"]

    server = "blast.ncbi.nlm.nih.gov"
    url = "/blast/Blast.cgi"
    try:#{
        conn = http.client.HTTPSConnection(server) # create a connection
        conn.request("POST", url, payload, headers) # send the request
        response = conn.getresponse() # get the response
    #}
    except OSError as oserr:#{
        print(get_work_time() + "\n - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
        print("Check your Internet connection.\a")
        print('\n' + '=/' * 20)
        print( repr(err) )
        platf_depend_exit(-2)
    #}

    response_text = str(response.read(), "utf-8") # get response text

    # /===== TEST =====\
    # with open("response.txt", 'w') as resp_file:
    #     resp_file.write(response_text + '\n')
    # \================/

    rid = re_search("RID = (.*)", response_text).group(1) # get Request ID
    rtoe = int(re_search("RTOE = (.*)", response_text).group(1)) # get time to wait provided by the NCBI server

    conn.close()

    return {"RID": rid, "RTOE": rtoe}
#}


def wait_for_align(rid, rtoe, attempt, attempt_all, filename):#{
    """
    Function waits for satisfaction of the request and retrieves response text.

    :param rid: Reques ID to wait for;
    :type rid: str;
    :param rtoe: time in seconds to wait provided by the NCBI server;
    :type rtoe: int;
    :param attempt: current attempt. This information is printed to console;
    :type attempt: int;
    :param attempt all: total number of attempts corresponding to current FASTA file.
        This information is printed to console;
    :type attempt_all: int;
    :param filename: basename of current FASTA file;
    :type filename: str;

    Returns XML text of type 'str' with BLAST response.
    """
    
    
    print("\n{} - Requested data for: '{}' ({}/{})".format(get_work_time(), filename, attempt, attempt_all))
    if rtoe > 0:#{ RTOE can be zero at the very beginning of continuation
        print("{} - The system estimated that the query with Request ID '{}' will be resolved in {} seconds ".format(get_work_time(), rid, rtoe))
        print("{} - Going to sleep for that period...".format(get_work_time()))
        # Server migth be wrong -- we will give it 3 extra seconds
        sleep(rtoe + 3)
        print("{} - Woke up. Checking request updates...".format(get_work_time()))
    #}

    server = "blast.ncbi.nlm.nih.gov"
    wait_url = "/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid
    there_are_hits = False

    while True:#{
        # sleep(60) ???
        error = True
        while error:#{
            try:
                conn = http.client.HTTPSConnection(server) # create connection
                conn.request("GET", wait_url) # ask for if there areresults
            except TimeoutError as err:
                print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
                print('\t' + repr(err))
                error = True
                sleep(30)
            except http.client.RemoteDisconnected as err:
                print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
                print('\t' + repr(err))
                error = True
                sleep(30)
            except ConnectionResetError as err:
                print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
                print('\t' + repr(err))
                error = True
                sleep(30)
            except FileNotFoundError as err:
                print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
                print('\t' + repr(err))
                error = True
                sleep(30)
            except http.client.CannotSendRequest as err:
                print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
                print('\t' + repr(err))
                error = True
                sleep(30)
            else:
                error = False # if no exception ocured
        #}

        response = conn.getresponse() # get the resonse
        resp_content = str(response.read(), "utf-8") # get response text
        conn.close()

        if re_search("Status=WAITING", resp_content) is not None:#{ if server asks to wait
            print("{} - Still searching for '{}'.".format(get_work_time(), filename), end="")
            print(" Going to sleep for another 60 seconds.")
            sleep(60)
            continue
        #}
        if re_search("Status=FAILED", resp_content) is not None:#{ if query failed
            print('\n' + get_work_time() + " - Query failed\a\n")
            response_text = """{} - Query for {} with Request ID {} failed.
    Contact NCBI or try to start it again.\n""".format(get_work_time(), filename, rid)
            return None
        #}
        if re_search("Status=UNKNOWN", resp_content) is not None:#{ if query expired
            print('\n' + get_work_time() + " - Query expired\a\n")
            respond_text = """{} - Query for {} with Request ID {} expired.
    Try to start it again\n""".format(get_work_time(), filename, rid)
            return "expired"
        #}
        if re_search("Status=READY", resp_content) is not None:#{ if results are ready
            there_are_hits = True
            print("\n{} - Result for query '{}' ({}/{}) is ready!".format(get_work_time(), filename, attempt, attempt_all))
            if re_search("ThereAreHits=yes", resp_content) is not None:#{ if there are hits
                print(get_work_time() + " - There are hits. Retrieving them.")
                for i in range(45, 0, -5):
                    print('-' * i)
                break
            #}
            else:#{ if there are no hits
                print(get_work_time() + " - There are no hits. It happens.\n")
                break
            #}
        #}
        # Execution should not reach here
        print('\n' + get_work_time() + " - Something unexpected happened. Contact the developer.\a\n")
        platf_depend_exit(1)
    #}

    # Retrieve XML result
    retrieve_xml_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&ALIGNMENTS=1&RID=" + rid
    conn = http.client.HTTPSConnection(server)
    conn.request("GET", retrieve_xml_url)
    response = conn.getresponse()

    respond_text = str(response.read(), "utf-8")
    conn.close()

    # # /==== TEST ====\
    # with open("align_text.xml", 'w') as alfile:
    #     alfile.write(respond_text + '\n')
    # # \==============/


    # Retrieve human-readable text and put it into result directory
    if there_are_hits:#{

        retrieve_text_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&DESCRIPTIONS=1&ALIGNMENTS=1&RID=" + rid
        conn = http.client.HTTPSConnection(server)
        conn.request("GET", retrieve_text_url)
        response = conn.getresponse()

        txt_hpath = "{}{}blast_result_{}.txt".format(filename[:filename.rfind(".fasta")], os.sep, attempt)
        # Write text result for a human to read
        with open(txt_hpath, 'w') as txt_file:
            txt_file.write(str(response.read(), "utf-8") + '\n')
        conn.close()
    #}

    return respond_text
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

            pident_ratio = round( float(pident) / int(align_len) * 100, 2)
            gaps_ratio = round( float(gaps) / int(align_len) * 100, 2)

            print("""\n\tHit!: '{}' -- '{}'
                Identity - {}/{} ({}%); Gaps - {}/{} ({}%)""".format(query_name, hit_taxa_name,
                    pident, align_len, pident_ratio, gaps, align_len, gaps_ratio))

            # Append new tsv line containing recently collected information
            result_tsv_lines.append( DELIM.join( [query_name, hit_taxa_name, hit_acc,
                align_len, pident, gaps, evalue] ))
        #}
        else:
            # If there is no hit for current sequence
            print("\t{} -- No significant similarity found".format(query_name))
            result_tsv_lines.append(DELIM.join( [query_name, "No significant similarity found"] ))
        print('=' * 30 + '\n')
    #}
    return result_tsv_lines
#}


def write_result(res_tsv_lines, tsv_res_path):#{
    """
    Function writes result of blasting to result tsv file and sorts
        recently blasted reads between corresponding FASTQ files.

    :param res_tsv_lines: tsv lines returned by 'parse_align_results_xml()' funciton;
    :type res_tsv_lines: list<str>;
    :param tsv_res_path: path to reslut tsv file;
    :type tsv_res_path: str;
    """


    # If there is no result tsv fil -- create it and write a head of the table.
    if not os.path.exists(tsv_res_path):#{
        with open(tsv_res_path, 'w') as tsv_res_file:#{
            tsv_res_file.write(DELIM.join( ["QUERY_ID", "HIT_NAME", "HIT_ACCESSION",
                "ALIGNMENET_LENGTH", "IDENTITY(%)", "GAPS(%)", "E-VALUE"] ) + '\n')
        #}
    #}
    # Write reslut tsv lines to this file
    with open(tsv_res_path, 'a') as tsv_res_file:#{
        for line in res_tsv_lines:
            tsv_res_file.write(line + '\n')
    #}

#}


def create_result_directory(fastq_path):#{
    """
    Function creates a result directory named according to how source FASTQ file is named.

    :param fastq_path: path to source fastq file;
    :type fastq_path: str;

    Returns 'str' path to the recently created result directory.
    """

    # dpath means "directory path"
    new_dpath = os.path.basename(fastq_path) # get rid of absolute path
    new_dpath = new_dpath[: new_dpath.rfind(".fastq")] # get rid if '.fastq' or 'fastq.gz'
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
#        "attmpt": saved_attempt (int),
#        "RID": saved_RID (str),
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. 'response' is a dic of the following structure:
#    {
#        "RID": Request ID (str),
#        "RTOE", time_to_wait_provided_by_the_NCBI_server (int)
#    }

#                   |===== Kernel loop =====|

fastq_list = get_fastq_list() # get source FASTQ files to process
print("\nChecking Internet connection...")
check_connection()
print("OK\n")
organisms = get_organisms() # get organisms for 'nt' database restriction
blast_algorithm = get_algorithm() # get BLAST algorithm


# Iterate through found source FASTQ files
for i, fastq_path in enumerate(fastq_list):#{
    
    # Create the result directory with the name of FASTQ file being processed:
    new_dpath = create_result_directory(fastq_path)

    # Convert fastq file to FASTA and get it's path and number of reads in it:
    curr_fasta = fastq2fasta(fastq_path, i, new_dpath)

    # "hname" means human readable name (e.i. without file path and extention)
    fasta_hname = os.path.basename(curr_fasta["fpath"]) # get rid of absolure path
    fasta_hname = fasta_hname[: fasta_hname.rfind(".fasta")] # get rid of file extention

    # Look around and ckeck if there are results of previous runs of this script
    # If 'look_around' is None -- there is no data from previous run
    previous_data = look_around(new_dpath, curr_fasta["fpath"],
        blast_algorithm)

    if previous_data is None:#{ # If there is no data from previous run
        num_done_reads = 0 # number of successfully processed reads
        saved_attempt = None # number of last successful attempt (there is no such stuff for de novo run)
        packet_size = get_packet_size(curr_fasta["nreads"]) # ask a user for packet size
        tsv_res_path = "{}_{}_result.tsv".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form result tsv file path
        tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form temporary file path
        with open(tmp_fpath, 'w') as tmp_file:
            tmp_file.write(str(packet_size)+ '\n')
    #}
    else:#{ # if there is data from previous run
        num_done_reads = previous_data["n_done_reads"] # get number of successfully processed reads
        saved_attempt = previous_data["attmpt"] # get number of last successful attempt
        packet_size = previous_data["pack_size"] # packet size sholud be the same as it was in previous run
        tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        tmp_fpath = previous_data["tmp_fpath"] # temporary file sholud be the same as during previous run
        saved_RID = previous_data["RID"] # having this RID we can try to get response for last request
        contin_rtoe = 0 # we will not sleep at the very beginning of continuation
    #}

    attempt_all = curr_fasta["nreads"] // packet_size # Calculate total number of packets sent from current FASTA file
    if curr_fasta["nreads"] % packet_size > 0: # And this is ceiling (in order not to import 'math')
        attempt_all += 1
    attempts_done = int( num_done_reads / packet_size ) # number of successfully processed reads

    with open(curr_fasta["fpath"], 'r') as fasta_file:#{

        # Go untill the last processed sequence
        for _ in range( int(num_done_reads * 2) ):#{
            fasta_file.readline()
        #}

        reads_left = curr_fasta["nreads"] - num_done_reads # number of reads left to precess
        attempts_left = attempt_all - attempts_done # number of packets left to send
        attempt = attempts_done+1 if attempts_done > 0 else 1 # current attempt

        # Iterate througth packets left to send
        for i in range(attempts_left):#{

            packet = get_packet(fasta_file, packet_size) # form the packet

            if packet["fasta"] is "":#{   Just in case
                print("Recent packet is empty")
                break
            #}

            print("\nGo to BLAST (" + blast_algorithm + ")!")
            print("Request number {} out of {}.".format(attempt, attempt_all))

            # # /=== Test ===\
            # with open("align_text.xml", 'r') as xml_file:
            #     xml_text = xml_file.read()
            # result_tsv_lines = parse_align_results_xml(xml_text,
            #                 packet["names"])
            # write_result(result_tsv_lines, tsv_res_path)
            # continue
            # # \============/

            send = True

            # If current packet has been already send, we can try to get it and not to send it again
            if attempt == saved_attempt and saved_RID is not None:#{

                align_xml_text = wait_for_align(saved_RID, contin_rtoe,
                    attempt, attempt_all, fasta_hname+".fasta") # get BLAST XML response

                # If request is not expired get he result and not send it again
                if align_xml_text != "expired":#{
                    send = False

                    result_tsv_lines = parse_align_results_xml(align_xml_text,
                        packet["names"]) # get result tsv lines

                    # Write the result to tsv and sort FASTQ sequences
                    write_result(result_tsv_lines, tsv_res_path)
                #}
            #}

            if send:#{
            
                request = configure_request(packet["fasta"], blast_algorithm, organisms) # get the request
                response = send_request(request) # send the request

                # Save temporary data
                with open(tmp_fpath, 'a') as tmpfile:
                    tmpfile.write("{}\t{}\n".format(attempt, response["RID"]))

                align_xml_text = wait_for_align(response["RID"], response["RTOE"],
                    attempt, attempt_all, fasta_hname+".fasta") # get BLAST XML response

                result_tsv_lines = parse_align_results_xml(align_xml_text,
                    packet["names"]) # get result tsv lines

                # Write the result to tsv and sort FASTQ sequences
                write_result(result_tsv_lines, tsv_res_path)
            #}

            attempt += 1
        #}
    #}

    # Remove temporary file
    if os.path.exists(tmp_fpath):
        os.unlink(tmp_fpath)
#}


print("\nTask completed successfully!")
platf_depend_exit(0)