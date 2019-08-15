#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# |===== Check python interpreter version =====|

from sys import version_info as verinf

if verinf.major < 3:#{
    print( "Your python interpreter version is " + "%d.%d" % (verinf.major, verinf.minor) )
    print("\tPlease, use Python 3!\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3
    raw_input("Press ENTER to exit:")
    exit(1)
#}

# |===== Function for dealing with time =====|

from time import time, strftime, gmtime, localtime, sleep
start_time = time()

print( strftime("\n%H:%M:%S", localtime(start_time)) + " - START WORKING\n")

def get_work_time():#{
    return strftime("%H:%M:%S", gmtime( time() - start_time ))
#}

import os
from re import search as re_search
from gzip import open as open_as_gzip # input files might be gzipped
from xml.etree import ElementTree

import http.client
import urllib.request
from urllib.error import HTTPError
import urllib.parse


# |===== Function that will ask for ENTER pressing on Windows =====|

from sys import platform

def platf_depend_exit(exit_code):#{
    if "win" in sys.platform:
        input("Press ENTER to exit:")
    exit(exit_code)
#}


# |===== Function for checking if 'https://blast.ncbi.nlm.nih.gov' is available =====|

def check_connection():#{

    try:#{

        ncbi_server = "https://blast.ncbi.nlm.nih.gov"
        status_code = urllib.request.urlopen(ncbi_server).getcode()

        # Just in case
        if status_code != 200:#{
            print(get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            print("Check your Internet connection.\a")
            print("Status code: {}".format(status_code))
            platf_depend_exit(-2)
        #}
        return
    #}
    except OSError as err:#{

        print(get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
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
print("Checking Internet connection...")
check_connection()
print("OK\n" + "~"*30 + '\n')



# |===== Question-funtions =====|

# Функция для опроса пользователя, желает ли он продолжить работу скрипта
def is_continued():#{

    continuation = None

    while continuation is None:#{
        continuation = input("""
Would you like to continue the previous run?
    1. Continue!
    2. Start from the beginning.

Enter the number (1 or 2):>> """)

    # Проверка, введено ли число. Если нет, то даётся очередная попытка
        try:#{
            continuation = int(continuation)
            if continuation != 1 and continuation != 2:#{ Проверка, введено ли число от 1 до 2.
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


# Запрос размера пакета для обработки fastq-файлов
def get_packet_size(num_reads):#{

    packet_size = None
    limit = num_reads if num_reads <= 1000 else 1000

    while packet_size is None:#{
        
        packet_size = input("""
Please, specify the number of sequences that should be sent to the NCBI server in one request.
Enter the number (from 1 to {}):>> """.format(limit))
        # Проверка, введено ли число. Если нет, то даётся очередная попытка
        try:#{
            packet_size = int(packet_size)
            if packet_size < 1 or packet_size > limit:#{
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

    sens = None

    while sens is None:#{
        sens = input("""
Please, specify the taxonomy level to classify by:
    1. Genus.
    2. Species.
    3. Strain.
Enter a number (1, 2 or 3):>> """)
        try:#{
            sens = int(sens)
            if sens < 1 or sens > 3:#{
                print("\n\tNot a valid number entered!\a\n" + '~'*20)
                sens + None
            #}
            else:#{
                print("You have chosen number "+ str(sens) + '\n')
                print('~' * 20 + '\n')

                if sens is 1:
                    sens = "genus"
                elif sens is 3:
                    sens = "strain"
                else:
                    sens = "species"   # this behaviour is some kind of default
                
            #}
        #}
        except ValueError:#{
            print("\nNot an integer NUMBER entered!\a\n" + '~'*20)
            algorithm = None
        #}
    #}
    return sens
#}

# Запрос выбора алгоритма BLAST для обработки fasta-файлов
def get_algorithm():#{

    reply = None
    blast_algorithm = "megaBlast" # default value

    while reply is None:#{
        reply = input("""
Please choose a BLAST algorithm:
    1. Highly similar sequences (megablast)
    2. Optimize for More dissimilar sequences (discontiguous megablast)
    3. Optimize for Somewhat similar sequences (blastn)

Enter the number (1 or 2 or 3):>> """)
        # Проверка, введено ли число. Если нет, то даётся очередная попытка
        try:#{
            reply = int(reply)
            if reply < 1 or reply > 3:#{ Проверка, введено ли число от 1 до 3.
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


def get_fastq_list():#{

    # Составляем список файлов для обработки - с расширением ".fastq"
    fastq_list = os.listdir('.')  # в переменную fastq_list записываем все названия файлов из рабочей директории
    # Затем оставляем только файлы, заканчиващиеся на ".fastq" или ".fastq.gz"
    is_fastq_or_fastqgz = lambda f: f.endswith(".fastq") | f.endswith(".fastq.gz")
    fastq_list = list(filter(is_fastq_or_fastqgz, fastq_list))

    if len(fastq_list) != 0:#{
        print("Following files have been found and will be processed:")
        for i, line in enumerate(fastq_list):#{
            print("\t{}. {}".format( str(i+1), fastq_list[i] ))
        #}
        print('-' * 40)
    #}
    else:#{
        print("No '.fastq' or '.fastq.gz' files have been found in current directory!")
        platf_depend_exit(1)  # завершаем
    #}
    return fastq_list
#}


# |===== Functionality for proper processing of gzipped files =====|

OPEN_FUNCS = (open, open_as_gzip)

# Assuming that people wouldn't mislead anybody by file renaming
is_gzipped = lambda file: True if file.endswith(".gz") else False

# Data from .fastq and .fastq.gz should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format .fastq line
    lambda line: line.decode("utf-8").strip()  # format .fastq.gz line
)



FASTQ_LINES_PER_READ = 4
FASTA_LINES_PER_READ = 2

def fastq2fasta(fastq_path, i):#{
    """
    :param i: order number of fastq_file
    """

    how_to_open = OPEN_FUNCS[is_gzipped(fastq_path)]
    fmt_func = FORMATTING_FUNCS[is_gzipped(fastq_path)]
    
    # dpath means "directory path"
    new_dpath = os.path.basename(fastq_path) # get rid of absolute path
    new_dpath = new_dpath[: new_dpath.rfind(".fastq")] # get rid if '.fastq' or 'fastq.gz'
    
    fasta_path = fastq_path.replace(".fastq", ".fasta")
    fasta_path = os.path.join(new_dpath, fasta_path)

    if not os.path.exists(new_dpath):#{
        os.makedirs(new_dpath)
    #}

    num_lines = 0
    if not os.path.exists(fasta_path):#{
        with how_to_open(fastq_path) as fastq_file, open(fasta_path, 'w') as fasta_file:#{

            counter = 1
            for line in fastq_file:#{  считываем файл построчно
                line = fmt_func(line)
                if counter <= 2:#{      записываем только первые две строки из четырёх
                    if line[0] == '@':
                        line = '>' + line[1:]  # заменяем собаку на знак "больше" - как и положено в fasta-формате
                    fasta_file.write(line + '\n')
                #}
                elif counter == 4:
                    counter = 0
                counter += 1
                num_lines += 1
            #}
        #}
        num_reads = int(num_lines / FASTQ_LINES_PER_READ)
    #}
    else:#{
        num_lines = sum(1 for line in open(fasta_path, 'r'))
        num_reads = int(num_lines / FASTA_LINES_PER_READ)
    #}

    print("{}. '{}' ({} read) --> FASTA\n".format(i+1, fastq_path, num_reads))

    return {"dpath": new_dpath, "fpath": fasta_path, "nreads": num_reads}
#}


def peek_around(new_dpath, fasta_path, blast_algorithm):#{
    

    # "hname" means human readable name (e.i. without file path and extention)
    fasta_hname = os.path.basename(fasta_path) # get rid of absolute path
    fasta_hname = fasta_hname[: fasta_hname.rfind(".fasta")] # get rid of '.fasta' extention

    tmp_fpath = "{}.{}_temp.txt".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)
    tsv_res_fpath = "{}.{}_result.tsv".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)

    # Отправляем полученный fasta-файл на сервер NCBI для BLASTn анализа
    print("========== file: '{}' ===========".format(fasta_path))
    continuation = False    # По умолчанию делаем значение "продолжения" на "de novo"
    # Проверяем, есть ли файлы вывода с предыдущего раза, предоставляем возможность продолжить с последнего успешного результата
    if os.path.exists(tmp_fpath):#{
        print(get_work_time() + " - The previous result file is found in the directory:")
        print("\t'{}'".format(tmp_fpath))
        continuation = is_continued()
    #}

    if continuation:#{   # Находим название последней проанализированной последовательности
        print("Let\'s try to continue...")
        if os.path.exists(tsv_res_fpath):#{
            with open(tsv_res_fpath, 'r') as res_file:#{
                lines = res_file.readlines()
                num_done_reads = len(lines)
                last_line = lines[-1]
                last_seq_id = last_line.split('\t')[0]
            #}
            print("Last successful attempt: " + last_seq_id)
        #}
        else:#{
            continuation = False                                                               # WHAT FOR ???
        #}

        try:
            _ = num_done_reads
        except NameError:
            num_done_reads = 0

        # Находим номер последней попытки и RID, а также размер пакета
        with open(tmp_fpath, 'r') as tmp_file:
            temp_lines = tmp_file.readlines()
        packet_size = int(temp_lines[0])
        attempt_save = int(temp_lines[-1].split('\t')[0])
        RID_save = temp_lines[-1].split('\t')[1].strip()

        return {
            "pack_size": packet_size,
            "attmpt": attempt_save,
            "RID": RID_save,
            "tsv_respath": tsv_res_fpath,
            "n_done_reads": num_done_reads,
            "tmp_fpath": tmp_fpath
        }

    else:#{ Удаляем файлы вывода с предыдущего запуска сценария, если запускаем поиск от начала файла
        if os.path.exists(tsv_res_fpath):#{
            print(get_work_time() + " - The following file is found in the current directory and will be deleted:")
            print("\t'{}'".format(tsv_res_fpath))
            os.unlink(tsv_res_fpath)
        #}
        if os.path.exists(tmp_fpath):#{
            print(get_work_time() + " - The following file is found in the current directory and will be deleted:")
            print("\t'{}'".format(tmp_fpath))
            os.unlink(tmp_fpath)
        #}

        return None
    #}
#}


def get_packet(fasta_file, packet_size):#{

    packet = ""
    names = list()

    for i in range(packet_size):#{

        seq_id = fasta_file.readline().strip().partition(' ')[0] # prune the seq id a little bit
        names.append(seq_id)
        # seq_id = fasta_file.readline().strip()
        seq = fasta_file.readline().strip()

        # 'seq' cannot be empty string
        if seq_id == "":  # if previous sequence was last in current file
            return packet.strip()  # remove the last '\n' character
        else:
            packet += "{}\n{}\n".format(seq_id, seq)
    #}

    return {"fasta": packet.strip(), "names": names}  # remove the last '\n' character
#}


def configure_request(packet, blast_algorithm):#{

    query = packet
    database = "nt"
    taxonomy = "Escherichia (taxid:561)"
    # taxonomy = "bacteria (taxid:2)" # Необязательные параметры для среза БД nt
    # taxonomy1 = "viruses (taxid:10239)" # Необязательные параметры для среза БД nt
    # taxonomy = "Pseudomonas (taxid:286)" # Необязательные параметры для среза БД nt
    taxonomy1 = ""
    # taxonomy1 = "Rhodococcus (taxid:1827)" # Необязательные параметры для среза БД nt
    num_org = '2' # Необязательные параметры для среза БД nt
    program = "blastn"

    if "megablast" in blast_algorithm.lower():
        megablast = "on"
    else:
        megablast = ""

    cmd = "PUT"
    payload = {
        "CMD" : cmd,
        "PROGRAM": program,
        "MEGABLAST" : megablast,
        "BLAST_PROGRAMS" : blast_algorithm,
        "DATABASE" : database,
        "EQ_MENU" : taxonomy,   # Необязательные параметры
        "EQ_MENU1" : taxonomy1, # Необязательные параметры
        "NUM_ORG" : num_org,    # Необязательные параметры
        "QUERY" : query,
        "HITLIST_SIZE": 1
    }
    payload = urllib.parse.urlencode(payload)
    headers = {
        "Content-Type" : "application/x-www-form-urlencoded"
    }

    return {"payload":payload, "headers": headers}
#}

def send_request(request):#{
    
    try:#{
        payload = request["payload"]
        headers = request["headers"]
    #}
    except KeyError as kerr:#{
        print(repr(kerr))
        print("KeyError in 'send_request' function\a")
        platf_depend_exit(1)
    #}

    server = "blast.ncbi.nlm.nih.gov" #HTTPS
    url = "/blast/Blast.cgi"
    try:#{
        conn = http.client.HTTPSConnection(server)
        conn.request("POST", url, payload, headers)
        response = conn.getresponse()
    #}
    except OSError as oserr:#{
        print(get_work_time() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
        print("Check your Internet connection.\a")
        print('\n' + '=/' * 20)
        print( repr(err) )
        platf_depend_exit(-2)
    #}

    response_text = str(response.read(), "utf-8")

    # with open("response.txt", 'w') as resp_file:
    #     resp_file.write(response_text + '\n')

    rid = re_search("RID = (.*)", response_text).group(1)
    rtoe = int(re_search("RTOE = (.*)", response_text).group(1))

    conn.close()

    return {"RID": rid, "RTOE": rtoe}
#}


def wait_for_align(rid, rtoe, attempt, attempt_all, filename):#{
    
    try:
        _ = rtoe                                                                    # WHAT FOR ???
    except NameError:
        print("\n{} - Looking for respond for: {} ({}/{})".format(get_work_time(), filename, attempt, attempt_all))
        print("{} - Checking query with Request ID {}.".format(get_work_time(), rid))
    else:
        print("\n{} - Requested data for: '{}' ({}/{})".format(get_work_time(), filename, attempt, attempt_all))
        print("{} - The system estimated that the query with Request ID '{}' will be resolved in {} seconds ".format(get_work_time(), rid, rtoe))
        print("{} - Going to sleep for that period...".format(get_work_time()))
        # Server migth be wrong -- we will give it 3 extra seconds
        sleep(rtoe + 3)
        print("{} - Woke up. Checking request updates...".format(get_work_time()))

    server = "blast.ncbi.nlm.nih.gov" #HTTPS
    wait_url = "/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid
    there_are_hits = False

    while True:#{
        sleep(60)
        error = True
        while error:#{
            try:
                conn = http.client.HTTPSConnection(server)
                conn.request("GET", wait_url)
                error = False
            except TimeoutError as err:
                print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
                print('\t' + repr(err))
                error = True
                sleep(30)
            # except socket.gaierror as err:
            #     print("{} - Unable to connect to the NCBI server. Let\'s try to connect in 30 seconds.".format(get_work_time()))
            #     print('\t' + repr(err))
            #     error = 0
            #     sleep(30)
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
        #}

        response = conn.getresponse()
        resp_content = str(response.read(), "utf-8")
        conn.close()

        if re_search("Status=WAITING", resp_content) is not None:#{
            print("{} - Still searching for '{}'.".format(get_work_time(), filename))
            print("\tGoing to sleep for another 60 seconds.")
            continue
        #}
        if re_search("Status=FAILED", resp_content) is not None:#{
            print(get_work_time() + " - Query failed\a\n")
            response_text = """{} - Query for {} with Request ID {} failed.
    Contact NCBI or try to start it again.\n""".format(get_work_time(), filename, rid)
            return None
        #}
        if re_search("Status=UNKNOWN", resp_content) is not None:#{
            print(get_work_time() + " - Query expired\a\n")
            respond_text = """{} - Query for {} with Request ID {} expired.
    Try to start it again\n""".format(get_work_time(), filename, rid)
            return "expired"
        #}
        if re_search("Status=READY", resp_content) is not None:#{
            there_are_hits = True
            print("{} - Result for query '{}' ({}/{}) is ready!".format(get_work_time(), filename, attempt, attempt_all))
            if re_search("ThereAreHits=yes", resp_content) is not None:#{
                print(get_work_time() + " - There are hits. Retrieving them.")
                for i in range(45, 0, -5):
                    print('-' * i)
                break
            #}
            else:#{
                print(get_work_time() + " - There are no hits. It happens.\n")
                break
            #}
        #}
        print(get_work_time() + " - Something unexpected happened. Contact the developer.\a\n")
        platf_depend_exit(1)
        
    #}

    # Retrieve XML result
    retrieve_xml_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&ALIGNMENTS=1&RID=" + rid # + Tabular, - DESCRIPTIONS
    conn = http.client.HTTPSConnection(server)
    conn.request("GET", retrieve_xml_url)
    response = conn.getresponse()

    respond_text = str(response.read(), "utf-8")
    conn.close()


    # /==== TEST ====\
    with open("align_text.xml", 'w') as alfile:
        alfile.write(respond_text + '\n')
    # \==============/

    if there_are_hits:#{
        # Retrieve human-readable text
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

#                      Genus          species (mb sp.)   strain id  anything after comma
hit_name_pattern = r"^[A-Z][a-zA-Z]+ [a-zA-Z]*(sp\.)? (strain )?.+, .+$"

def format_taxonomy_name(hit_name, sens):#{

    # If structure of hit name is unforseen
    if re_search(hit_name_pattern, hit_name.strip()) is None:#{
        print("\tAttention!")
        print("Hit name '{}' has structure, that is unforseen by the developer.".format(hit_name))
        print("This name migth be formatted incorrectly.")
        print("Full name ({}) will be used.".format(taxa_name))
        print("Contact the develooper -- send this name to him.")

        return taxa_name.replace(' ', '_')   # return full name
    #}

    taxa_name = hit_name.partition(',')[0]
    taxa_splitnames = taxa_name.strip().split(' ')
    if sens == "genus":#{
        return taxa_name[0] # return genus
    #}
    elif sens == "species":#{
        if taxa_splitnames[1] == "sp.":#{ if there is no species specified
            return taxa_name.replace(' ', '_')   # return full name
        #}
        else:#{
            return '_'.join(taxa_splitnames[0], taxa_splitnames[1]) # return genus and species
        #}
    #}
    elif sense == "strain":#{
        return taxa_name.replace(' ', '_')   # return full name
    #}

    # Execution should not reach here
    raise Exception("Taxonomy nama formatting error!")
#}


def parse_align_results_xml(xml_text, seq_names, sens):#{

    result_tsv_lines = list()

    if "Bad Gateway" in xml_text:#{
        print('\n' + '=' * 45)
        print(get_work_time() + " - ERROR! Bad Gateway! Data from last pocket has lost.")
        print("It would be better if you restart the script.")
        print("Here are names of lost queries:")
        for i, name in enumerate(seq_names):#{
            print("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + "\t- Query has been lost: ERROR, Bad Gateway\n")
        #}
        input("Press ENTER to continue...")

        return result_tsv_lines
    #}

    if "to start it again" in xml_text:#{
        print(get_work_time() + "BLAST ERROR!")

        print("Here are names of lost queries:")
        for i, name in enumerate(seq_names):#{
            print("{}. '{}'".format(i+1, name))
            result_tsv_lines.append(name + "\t- Query has been lost: BLAST ERROR\n")
        #}

        input("Press ENTER to continue...")

        return result_tsv_lines
    #}

    # /=== Parse BLAST XML response ===/
    root = ElementTree.fromstring(xml_text)

    for iter_elem, iter_hit in zip(root.iter("Iteration"), root.iter("Iteration_hits")):#{
    
        query_name = iter_elem.find("Iteration_query-def").text

        hit = iter_hit.find("Hit")
        if hit is not None:#{

            hit_name = hit.find("Hit_def").text
            hit_taxa_name = format_taxonomy_name(hit_name, sens)
            print("'{}' -- {}".format(query_name, hit_taxa_name))

            hit_acc = hit.find("Hit_accession").text

            # Find the first HSP (we need only the first one)
            hsp = next(hit.find("Hit_hsps").iter("Hsp"))

            align_len = hsp.find("Hsp_align-len").text.strip()

            pident = float( hsp.find("Hsp_identity").text )
            pident = round( pident / align_len, 2 )

            gaps = float( hsp.find("Hsp_gaps").text )
            gaps = round( gaps / align_len, 2 )

            evalue = hsp.find("Hsp_evalue").text

            result_tsv_lines.append( '\t'.join(query_name, hit_taxa_name, hit_acc,
                align_len, pident, gaps, evalue, '\n') )
        #}
        else:
            result_tsv_lines.append('\t'.join(query_name, "No significant similarity found\n"))
    #}
    return result_tsv_lines
#}


def write_result(res_tsv_lines, res_dpath, source_fastq_path, tsv_res_path):#{

    with open(tsv_res_path, 'a') as tsv_res_file:#{
        for line in result_tsv_lines:
            tsv_res_file.write(line)
    #}

    how_to_open = OPEN_FUNCS[is_gzipped(fastq_path)]
    fmt_func = FORMATTING_FUNCS[is_gzipped(fastq_path)]

    for line in result_tsv_lines:#{

        seq_id = line.split['\t'][0]
        classfied_name = line.split['\t'][1]

        with how_to_open(source_fastq_path, 'r') as source_file:#{
            while True:#{
                fqline = fmt_func(source_file.readline())

                if fqline is "":#{
                    raise Exception("ERROR! Sent query sequence not found in source FASTQ file")
                #}

                if fqline.startswith('>'+seq_id):#{
                    id_line = fqline
                    seq = fmt_func(source_file.readline())
                    opt_id = fmt_func(source_file.readline())
                    qual_line = fmt_func(source_file.readline())
                    break
                #}
            #}
        #}

        classified_path = "{}{}{}_{}.fastq".format(res_dpath. os.sep, os.path.basename(fastq_path),
            classfied_name)

        with open(classified_path, 'a') as classif_file:#{
            for fqline in (id_line, seq, opt_id, qual_line):#{
                classif_file.write(fqline + '\n')
            #}
        #}
    #}
#}




#                |===== Proceed =====|
# =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=

fastq_list = get_fastq_list()
blast_algorithm = get_algorithm()
sensibility = get_classif_sensibility()


# /=== Comments to the kernel loop ===/

# 1. 'curr_fasta' is a dict of the following structure:
#    {
#       "dpath": path_to_directory_with_result_files,
#       "fpath": path_to_fasta_file,
#       "nreads": number_of_reads_in_this_fasta_file
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. 'previous_data' is a dict of the following structure:
#    {
#       "pack_size": packet_size,
#       "attmpt": saved_attempt,
#       "RID": saved_RID,
#       "tsv_respath": path_to_tsv_file_from_previous_run,
#       "n_done_reads": number_of_successfull_requests_from_currenrt_fasta_file,
#       "tmp_fpath": path to pemporary file
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# |===== Kernel loop =====|
for i, fastq_path in enumerate(fastq_list):#{

    # Convert fastq file to fasta and get it's path and number of reads in it:
    curr_fasta = fastq2fasta(fastq_path, i)

    # "hname" means human readable name (e.i. without file path and extention)
    fasta_hname = os.path.basename(curr_fasta["fpath"]) # get rid of absolure path
    fasta_hname = fasta_hname[: fasta_hname.rfind(".fasta")] # get rid of file extention

    # Peek around and ckeck if there are results of previous runs of this script
    # If 'peek_around' is None -- there is no data from previous run
    previous_data = peek_around(curr_fasta["dpath"], curr_fasta["fpath"], blast_algorithm)
    if previous_data is not None:#{
        num_done_reads = previous_data["n_done_reads"]
        packet_size = previous_data["pack_size"]
        saved_attempt = previous_data["attmpt"]
        saved_RID = previous_data["RID"]
        tsv_res_path = previous_data["tsv_respath"]
        tmp_fpath = previous_data["tmp_fpath"]
    #}
    else:#{
        num_done_reads = 0
        saved_attempt = None
        packet_size = get_packet_size(curr_fasta["nreads"])
        tmp_fpath = "{}.{}_temp.txt".format(os.path.join(curr_fasta["dpath"],
            fasta_hname), blast_algorithm)
        tsv_res_path = "{}.{}_result.tsv".format(os.path.join(curr_fasta["dpath"],
            fasta_hname), blast_algorithm)
        # Создаём временный файл для хранения значения размера пакета и RID-запросов
        with open(tmp_fpath, 'w') as tmp_file:
            tmp_file.write(str(packet_size)+ '\n')
    #}

    attempt_all = curr_fasta["nreads"] // packet_size # Подсчитываем число отправок пакетов на сервер NCBI с одного файла
    if curr_fasta["nreads"] % packet_size > 0: # А это округление в большую сторону, чтобы не использовать math
        attempt_all += 1
    attempts_done = int( num_done_reads / packet_size )

    with open(curr_fasta["fpath"], 'r') as fasta_file:#{

        # Добираемся до последней обработанной последовательности
        for _ in range( int(num_done_reads * 2) ):#{
            fasta_file.readline()
        #}

        reads_left = curr_fasta["nreads"] - num_done_reads
        attempts_left = attempt_all - attempts_done
        attempt = attempts_done if attempts_done > 0 else 1

        for i in range(attempts_left):#{

            packet = get_packet(fasta_file, packet_size)

            if packet is "":#{   Just in case
                print("Well done!")
                break
            #}

            print("\nGo to BLAST (" + blast_algorithm + ")!")
            print("Request number {} out of {}.".format(attempt, attempt_all))

            send = True

            if attempt == saved_attempt:#{

                align_xml_text = wait_for_align(response["RID"], response["RTOE"],
                    attempt, attempt_all, os.path.basename(curr_fasta["fpath"]))
                if align_xml_text != "expired":#{

                    send = False

                    result_tsv_lines = parse_align_results_xml(align_xml_text,
                        packet["names"], sensibility)

                    write_result(result_tsv_lines, curr_fasta["dpath"], fastq_path, tsv_res_path)
                #}
            #}

            if send:#{
            
                request = configure_request(packet["fasta"], blast_algorithm)
                response = send_request(request)

                # Save temporary data
                with open(tmp_fpath, 'a') as tmpfile:
                    tmpfile.write("{}\t{}\n".format(attempt, response["RID"]))

                align_xml_text = wait_for_align(response["RID"], response["RTOE"],
                    attempt, attempt_all, os.path.basename(curr_fasta["fpath"]))

                result_tsv_lines = parse_align_results_xml(align_xml_text,
                    packet["names"], sensibility)

                write_result(result_tsv_lines, curr_fasta["dpath"], fastq_path, tsv_res_path)
            #}

            attempt += 1
        #}
    #}

    if os.path.exists(tmp_fpath):#{
        os.unlink(tmp_fpath)
    #}
#}


print("\nOK")