#!/usr/bin/env python3

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
import zipfile # for getting taxid information for 'nt' database restriction

import http.client
import urllib.request
from urllib.error import HTTPError
import urllib.parse


# |===== Function that asks to press ENTER press on Windows =====|

from sys import platform

def platf_depend_exit(exit_code):#{
    """
    A function that asks to press ENTER press on Windows
        and exits.

    :type exit_code: int;
    """
    if "win" in platform:
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
            print(get_work_time() + "\n - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            print("Check your Internet connection.\a")
            print("Status code: {}".format(status_code))
            platf_depend_exit(-2)
        #}
        return
    #}
    except OSError as err:#{

        print(get_work_time() + "\n - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
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
print("\nChecking Internet connection...", end="")
check_connection()
print("OK\n" + "~"*30 + '\n')


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

    # Проверка, введено ли целое число. Если нет, то даётся очередная попытка
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


def get_classif_sensibility():
    """
    Function asks a user about the sensibility of classification.
    It means that, for example, if the user decides 'Genus' variant,
        sequences will be sorted by genus name -- species and strain names
        won't be taken into consideration.

    :return: one of the following strings: "genus", "species", "strain";
    :return type: str;
    """
    sens = None

    while sens is None:#{
        sens = input("""
Please, specify the taxonomy level to classify by:
    1. Genus.
    2. Species.
    3. Strain.
Enter a number (1, 2 or 3):>> """)
        # Check if entered value is integer number. If no, give another attempt.
        try:#{
            sens = int(sens)
            if sens < 1 or sens > 3:#{ Check if input number is in [1, 3]
                print("\n\tNot a valid number entered!\a\n" + '~'*20)
                sens + None
            #}
            else:#{
                print("You have chosen number "+ str(sens) + '\n')
                print('~' * 20 + '\n')

                if sens is 1:
                    sens = "genus"
                elif sens is 2:
                    sens = "species"
                else:
                    sens = "strain"
            #}
        #}
        except ValueError:#{
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            algorithm = None
        #}
    #}
    return sens
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
            algorithm = None
    #}
    return blast_algorithm
#}

taxid_zip_path = "new_taxdump.zip"

def get_organisms():#{
    
    global taxid_zip_path

    # NCBI keeps this information in zipfile of 99Mb.
    # There fore we will download it only once during the run of script,
    #    save to current directory, work with fetched zip file
    #    and remove it after successful end.

    # Downloading of taxid zip file
    if not os.path.exists(taxid_zip_path):#{

        import threading
        taxid_zip_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"

        def download_waiter():#{

            while not os.path.exists(taxid_zip_path):
                sleep(1)
            print()
            fsize = round(os.path.getsize(taxid_zip_path) / (1024**2), 1) # get megabytes
            expect_size, ln_len = 98.0, 50 # file no ftp site is of 98,5 MB

            while fsize < expect_size:#{
                fsize = round(os.path.getsize(taxid_zip_path) / (1024**2), 1)
                eqs = int( (fsize/expect_size) / (100/ln_len) * 100 )
                spcs = ln_len - eqs
                arr = '>' if eqs < ln_len else ""
                print('\r' + '[' + '='*eqs + arr + ' '*spcs + ']' + "({}/{}) MB ".format(fsize, expect_size), end="")
            #}
            print()
        #}

        waiter = threading.Thread(target=download_waiter)
        waiter.start()

        print("There is to taxid file yet")
        print("Downloading '{}'...".format(taxid_zip_url))
        
        urllib.request.urlretrieve(taxid_zip_url, taxid_zip_path)

        waiter.join()
        print("Downloading is completed\n\n")

        if not zipfile.is_zipfile(taxid_zip_path):#{
            print("ERROR: recently downloaded zip file '{}' is not a valid zipfile".format(taxid_zip_path))
            platf_depend_exit(1)
        #}
    #}

    # ~~~~~~~~~~~
    
    # Only 2 'nt' database restrictions for now
    max_org = 2
    organisms = list()
    
    error = True
    while error:#{
        reply = input("""~~~~~~~~~~~~~~~~~~~~~~~~~~~
Please specify organisms for 'nt' database restriction.
Enter taxid OR organism name after prompt '>>' symbol. For example:
\t>> 551
or
\t>> Erwinia

You can enter 2 taxids OR 2 organism names divided by comma as follows:
\t>> 551,1833
or
\t>> Erwinia,Rhodococcus erythropolis
or even
\t>> 551,Rhodococcus erythropolis
(551 is taxid of Erwinia, 1833 -- of Rhodocuccus erythropolis)

Default settings are: 'bacteria (taxid:2)' and 'viruses (taxid:10239)'.
To use default settings just press ENTER.

Enter organism name or taxid:>> """)

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
            org_list = reply.split(',')
            org_list = list( map(str.strip, org_list) )
            if len(org_list) > max_org:#{
                print("\n\t!!! -ERROR: You can enter only 1 or 2 organisms")
                print("This constraint is meant to be fixed\n")
                sleep(3)
                continue
            #}
        #}
        except Exception:#{
            print("Error: invalid input\n"+"~"*20+"\n")
            continue
        #}

        taxid_col = 0
        name_col = 1

        # Search for input entry in "new_taxdump.zip"
        print("Validating organisms...")
        for org in org_list:#{

            # If entrance can be interpreted as integer number -- probably, taxid was entered.
            # If it isn't -- probably, organism name was entered.
            try:
                _ = int(org)
            except ValueError:
                col_to_search = name_col # Search organism name (2-nd column in file)
            else:
                col_to_search = taxid_col # Search for taxid (1-st column in file)

            tax_zip =  zipfile.ZipFile(taxid_zip_path)
            names_dmp = tax_zip.open("names.dmp")

            org_found = False
            while not org_found:#{

                line = names_dmp.readline().decode("utf-8")
                if line == "":#{
                    print("\n\t!!! - ERROR: Organism(taxid) '{}' is not found\n\a".format(org))
                    sleep(3)
                    break
                #}
                check_obj = line.split('|')[col_to_search].strip()

                if check_obj.upper() == org.upper():#{
                    name = line.split('|')[name_col].strip()
                    taxid = line.split('|')[taxid_col].strip()

                    # It is better to choose scientific name
                    if col_to_search == 0:#{
                        tmp_taxid = taxid
                        # Search whrough names with this taxid (names with the same taxid are placed together in file)
                        while tmp_taxid == taxid:#{

                            line = names_dmp.readline().decode("utf-8")
                            if line == "":
                                break

                            tmp_taxid = line.split('|')[taxid_col].strip()
                            tmp_name = line.split('|')[name_col].strip()
                            comment = line.split('|')[3].strip()

                            if comment == "scientific name": # Find scientific name
                                name = tmp_name
                            #}
                    #}

                    format_tax = "{} (taxid:{})".format(name, taxid)
                    organisms.append(format_tax)
                    org_found = True
                #}
            #}
            tax_zip.close()
            names_dmp.close()
            if not org_found:
                organisms = list() # reset organisms list
                break
        #}
        error = True if not org_found else False
    #}
    print("\nYou have chosen following organisms:")
    for i, org in enumerate(organisms):#{
        print("\t{}. {}".format(i+1, org))
    #}
    print('~'*30 + '\n')
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
        print()
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
        num_reads = int(num_lines / FASTA_LINES_PER_READ) # get number of reads
    #}

    print("\n{}. '{}' ({} reads) --> FASTA\n".format(i+1, fastq_path, num_reads))

    return {"fpath": fasta_path, "nreads": num_reads}
#}


def remove_files_verbosely(*files):#{
    print('~' * 40)
    for file in files:#{
        try:#{
            print('\n' + get_work_time() + " - Following file will be removed:")
            print("\t'{}'".format(file))
            os.unlink(file)
        #}
        except:#{
            # Anything (and not only strings) can be passed to the function
            print("\nFile '{}' cannot be removed".format( str(file)) )
        #}
    #}
    print('~' * 40 + '\n')
#}


def peek_around(new_dpath, fasta_path, blast_algorithm):#{
    """
    Function peeks around in order to ckeck if there are results from previous runs of this script.

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
    tmp_fpath = "{}.{}_temp.txt".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)
    # Form path to result file
    tsv_res_fpath = "{}.{}_result.tsv".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)

    num_done_reads = None # variable to keep number of succeffdully processed reads

    print("\n |========== file: '{}' ===========|".format(fasta_path))
    continuation = False    # default value
    # Check if there are results from previous run.
    if os.path.exists(tsv_res_fpath):#{
        print('\n' + get_work_time() + " - The previous result file is found in the directory:")
        print("\t'{}'".format(tmp_fpath))
        continuation = is_continued() # Allow politely to continue from last successful attempt.
        if not continuation:
            remove_files_verbosely(tsv_res_fpath, tmp_fpath)
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
                    remove_files_verbosely(tsv_res_fpath, tmp_fpath)
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

    # /==== TEST ====\
    with open("align_text.xml", 'w') as alfile:
        alfile.write(respond_text + '\n')
    # \==============/


    # Retrieve human-readable text
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

#                      Genus    species                   strain name and anything after it
hit_name_pattern = r"^[A-Z][a-z]+ [a-z]*(sp\.)?(phage)? (strain )?.+$"

def format_taxonomy_name(hit_name, sens):#{
    """
    Function formats taxonomy name according to chosen sensibiliry of classification.

    :param hit_name: full_fit_name_of_the_subject_sequence;
    :type hit_name: str;
    :param sens: sensibility returned by 'get_classif_sensibility()' function.
        It's value can be one of the following strings: "genus", "sprcies", "strain";
    :type sens: str;

    Returns formatted hit name of 'str' type;
    """

    # If structure of hit name is strange
    if re_search(hit_name_pattern, hit_name.strip()) is None:#{
        print("\n\tAttention!")
        print("Hit name '{}' has structure that hasn't been anticipated by the developer.".format(hit_name))
        print("This name migth be formatted incorrectly.")
        print("Therefore, full name ({}) will be used.".format(taxa_name))
        print("Contact the develooper -- send this name to him.")

        return taxa_name.replace(' ', '_')   # return full name
    #}

    taxa_name = hit_name.partition(',')[0]
    taxa_splitnames = taxa_name.strip().split(' ')

    # If hit is a phage sequence
    if taxa_splitnames[1] == "phage":#{
        # Assumming that the man who sortes by genus or species isn't interested in phage strain name
        if sens == "genus" or sens == "species":
            return '_'.join( [taxa_splitnames[0], taxa_splitnames[1]] ) # return "<Host_name> phage"
        else:
            return taxa_name.replace(' ', '_')   # return full name if we sort by strain
    #}

    if sens == "genus":#{
        return taxa_splitnames[0] # return genus
    #}
    elif sens == "species":#{
        if taxa_splitnames[1] == "sp.":#{ if species is not specified
            return taxa_name.replace(' ', '_')   # return full name
        #}
        else:#{
            return '_'.join( [taxa_splitnames[0], taxa_splitnames[1]] ) # return genus and species
        #}
    #}
    elif sens == "strain":#{
        return taxa_name.replace(' ', '_')   # return full name
    #}

    # Execution should not reach here
    raise Exception("Taxonomy name formatting error!")
#}


def parse_align_results_xml(xml_text, seq_names, sens):#{
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
            # Format hit name accoding to classification sensitivity
            hit_taxa_name = format_taxonomy_name(hit_name, sens)

            hit_acc = hit.find("Hit_accession").text # get hit accession

            # Find the first HSP (we need only the first one)
            hsp = next(hit.find("Hit_hsps").iter("Hsp"))

            align_len = hsp.find("Hsp_align-len").text.strip()
            align_len = int(align_len) # we need 'align_len' as int for computations by now

            pident = float( hsp.find("Hsp_identity").text ) # get number of matched nucleotides
            pident = round( (pident / align_len) * 100, 2 ) # convert to percents and round it
            pident = str(pident) # we do not need it sa float any more

            gaps = float( hsp.find("Hsp_gaps").text ) # get number of gaps
            gaps = round( (gaps / align_len) * 100, 2 ) # convert to percents and round it
            gaps = str(gaps) # we do not need it sa float any more

            align_len = str(align_len) # we need 'align_len' as intany more

            evalue = hsp.find("Hsp_evalue").text # get e-value

            print("\n\tHit!: '{}' -- '{}' with e-value {}".format(query_name, hit_taxa_name, evalue))

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


def write_result(res_tsv_lines, res_dpath, source_fastq_path, tsv_res_path):#{
    """
    Function writes result of blasting to result tsv file and sorts
        recently blasted reads between corresponding FASTQ files.

    :param res_tsv_lines: tsv lines returned by 'parse_align_results_xml()' funciton;
    :type res_tsv_lines: list<str>;
    :param res_dpath: path to result directory;
    :type res_dpath: str;
    :param source_fastq_path: path to source FASTQ file;
    :type source_fastq_path: str;
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
        for line in result_tsv_lines:
            tsv_res_file.write(line + '\n')
    #}

    # Get ready for reading from '.fastq.gz'
    how_to_open = OPEN_FUNCS[ is_gzipped(fastq_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq_path) ]

    # Iretare through result tsv lines
    for line in result_tsv_lines:#{

        line = line.strip()

        seq_id = line.split(DELIM)[0] # seq id is the 1-st element of result tsv line

        # Find FASTQ record with recently got 'seq_id' in source FASTQ file
        with how_to_open(source_fastq_path, 'r') as source_file:#{
            while True:#{
                fqline = fmt_func(source_file.readline())

                # If there if no such read in source FASTQ file.
                # If it will occure, it is the developer mistake (not user's one).
                # Therefore -- raise an exception.
                if fqline is "":#{
                    raise Exception("ERROR! Sent query sequence not found in source FASTQ file")
                #}

                # If proper FASTQ record is found
                if fqline.startswith('@'+seq_id):#{
                    id_line = fqline # sequence id is recently read line
                    seq = fmt_func(source_file.readline()) # get sequence
                    opt_id = fmt_func(source_file.readline()) # get the 3-rd line
                    qual_line = fmt_func(source_file.readline()) # get qulity line
                    break # stop searching
                #}
            #}
        #}

        if "ERROR" in line:#{
            classified_name = "ERROR" # Erroneous reads will be written to separate file
        #}
        # Reads that perform no significant similarity will be written to separate file
        elif "No significant similarity found" in line:#{
            classified_name = "unknown"
        #}
        else:#{
            classified_name = line.split(DELIM)[1] # formatted hit name if the second element if tsv line
        #}

        # Form a result file path
        classified_path = "{}{}{}.{}.fastq".format(res_dpath, os.sep, os.path.basename(fastq_path),
            classified_name)
        # Colon in file name can provoke some problems on Windows
        classified_path = classified_path.replace(':', '')

        # Write FASTQ record to appropriate file
        with open(classified_path, 'a') as classif_file:#{
            for fqline in (id_line, seq, opt_id, qual_line):#{
                classif_file.write(fqline + '\n')
            #}
        #}
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
organisms = get_organisms() # get organisms for 'nt' database restriction
sensibility = get_classif_sensibility() # get sorting sensibility
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

    # Peek around and ckeck if there are results of previous runs of this script
    # If 'peek_around' is None -- there is no data from previous run
    previous_data = peek_around(new_dpath, curr_fasta["fpath"],
        blast_algorithm)

    if previous_data is None:#{ # If there is no data from previous run
        num_done_reads = 0 # number of successfully processed reads
        saved_attempt = None # number of last successful attempt (there is no such stuff for de novo run)
        packet_size = get_packet_size(curr_fasta["nreads"]) # ask a user for packet size
        tsv_res_path = "{}.{}_result.tsv".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form result tsv file path
        tmp_fpath = "{}.{}_temp.txt".format(os.path.join(new_dpath,
            fasta_hname), blast_algorithm) # form temporary file path
        # Создаём временный файл для хранения значения размера пакета и RID-запросов
        with open(tmp_fpath, 'w') as tmp_file:
            tmp_file.write(str(packet_size)+ '\n')
    #}
    else:#{ # if there is data from previous run
        num_done_reads = previous_data["n_done_reads"] # get number of successfully processed reads
        saved_attempt = previous_data["attmpt"] # get number of last successful attempt
        packet_size = previous_data["pack_size"] # packet size sholud be the same as during previous run
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
            #                 packet["names"], sensibility)
            # write_result(result_tsv_lines, new_dpath, fastq_path, tsv_res_path)
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
                        packet["names"], sensibility) # get result tsv lines

                    # Write the result to tsv and sort FASTQ sequences
                    write_result(result_tsv_lines, new_dpath, fastq_path, tsv_res_path)
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
                    packet["names"], sensibility) # get result tsv lines

                # Write the result to tsv and sort FASTQ sequences
                write_result(result_tsv_lines, new_dpath, fastq_path, tsv_res_path)
            #}

            attempt += 1
        #}
    #}

    # Remove temporary file
    if os.path.exists(tmp_fpath):
        os.unlink(tmp_fpath)
#}
# os.unlink(taxid_zip_path) # remove taxid-name mapping file


print("\nTask completed successfully!")
platf_depend_exit(0)