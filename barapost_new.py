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

from time import time, strftime, gmtime, localtime
start_time = time()

print( strftime("%H:%M:%S", localtime(start_time)) + " - START WORKING\n")

def get_work_time():#{
    return strftime("%H:%M:%S", gmtime( time() - start_time ))
#}

import os
# import re
from gzip import open as open_as_gzip # input files might be gzipped

import http.client
import urllib.request
from urllib.error import HTTPError
# import urllib.parse


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
# check_connection()
# print("OK\n" + "~"*30 + '\n')



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
            print("\nNot an integer NUMBER has been entered!\a\n" + '~'*20)
            continuation = None
        #}
    #}
    return(True if continuation == 1 else False)
#}


# Запрос размера пакета для обработки fastq-файлов
def get_packet_size():#{

    packet_size = None

    while packet_size is None:#{
        packet_size = input("""
Please, specify the number of sequences that should be sent to the NCBI server in one request.
Enter the number (from 1 to 1000):>> """)
        # Проверка, введено ли число. Если нет, то даётся очередная попытка
        try:#{
            packet_size = int(packet_size)
            if packet_size < 1 or packet_size > 1000:#{
                print("\n\tNot a VALID number entered!\a\n" + '~'*20)
                packet_size = None
            #}
            else:#{
                print("You have chosen number " + str(packet_size) + '\n')
                print('~' * 20 + '\n')
            #}
        #}
        except ValueError:#{
            print("\nNot an integer NUMBER has been entered!\a\n" + '~'*20)
            packet_size = None
        #}
    #}
    return(packet_size)
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
            print("\nNot an integer NUMBER has been entered!\a\n" + '~'*20)
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
    csv_res_fpath = "{}.{}_result.csv".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)

    # Отправляем полученный fasta-файл на сервер NCBI для BLASTn анализа
    print("========== file: '{}' ===========".format(fasta_path))
    continuation = False    # По умолчанию делаем значение "продолжения" на "de novo"
    # Проверяем, есть ли файлы вывода с предыдущего раза, предоставляем возможность продолжить с последнего успешного результата
    if os.path.exists(csv_res_fpath) and os.path.exists(tmp_fpath):#{
        print(get_work_time() + " - The previous result file is found in the directory:")
        print("\t'{}'".format(csv_res_fpath))
        continuation = is_continued()
    #}

    if continuation:#{   # Находим название последней проанализированной последовательности
        print("Let\'s try to continue...")
        if os.path.exists(csv_res_fpath):#{
            with open(csv_res_fpath, 'r') as res_file:#{
                lines = res_file.readlines()
                num_done_reads = len(lines)
                last_line = lines[-1]
                last_seq_id = last_line.split('\t')[0]
            #}
            print("Last successful attempt: " + last_seq_id)
        #}
        else:#{
            continuation = None                                                               # WHAT FOR ???
        #}

        # Находим номер последней попытки и RID, а также размер пакета
        with open(tmp_fpath, 'r') as tmp_file:
            temp_lines = tmp_file.readlines()
        packet_size = int(temp_lines[0])
        # attempt_save = int(temp_lines[-1].split('\t')[0])                                          # FOR NOW ????
        # RID_save = temp_lines[-1].split('\t')[1].strip()                                           # FOR NOW ????

        return {
        "pack_size": packet_size,
        # "attmpt": attempt_save,
        # "RID": RID_save,
        "csv_respath": csv_res_fpath,
        "n_done_reads": num_done_reads
        }

    else:#{ Удаляем файлы вывода с предыдущего запуска сценария, если запускаем поиск от начала файла
        if os.path.exists(csv_res_fpath):#{
            print(get_work_time() + " - The following file is found in the current directory and will be deleted:")
            print("\t'{}'".format(csv_res_fpath))
            os.unlink(csv_res_fpath)
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

    for i in range(packet_size):#{

        # seq_id = fasta_file.readline().strip().partition(' ')[0] # prune the seq id a little bit
        seq_id = fasta_file.readline().strip()
        seq = fasta_file.readline().strip()

        # 'seq' cannot be empty string
        if seq_id == "":  # if previous sequence was last in current file
            return packet.strip()  # remove the last '\n' character
        else:
            packet += "{}\n{}\n".format(seq_id, seq)
    #}

    return packet.strip()  # remove the last '\n' character
#}




#                |===== Proceed =====|
# =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=

fastq_list = get_fastq_list()
blast_algorithm = get_algorithm()



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
#       "csv_respath": path_to_csv_file_from_previous_run,
#       "n_done_reads": number_of_successfull_requests_from_currenrt_fasta_file
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
        packet_size = previous_data["pack_size"]
        # saved_attempt = previous_data["attmpt"]
        # saved_RID = previous_data["RID"]
        saved_csv_res_path = previous_data["csv_respath"]
        num_done_reads = previous_data["n_done_reads"]
    #}
    else:#{
        num_done_reads = 0
        packet_size = get_packet_size()
        tmp_fpath = "{}.{}_temp.txt".format(os.path.join(curr_fasta["dpath"],
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
        attempt = attempts_done


        print("Reads_at_all = {}".format(curr_fasta["nreads"]))
        print("Packet size = {}".format(packet_size))
        print("Attempt_all = {}".format(attempt_all))
        print("Attempts left = {}".format(attempts_left))
        if os.path.exists("test.fasta"):
            os.unlink("test.fasta")

        for i in range(attempts_left):#{

            print("Attempts done = {}".format(attempt))

            packet = get_packet(fasta_file, packet_size)

            if packet is "":#{  # Just in case
                print("Well done!")
                break
            #}

            print("\nGo to BLAST (" + blast_algorithm + ")!")
            print("Request number {} out of {}.".format(attempt+1, attempt_all))
            
            with open("test.fasta", 'a') as testfile:
                testfile.write(packet + '\n')
            # blastn_report = blast_the_packet(packet)


            attempt += 1
        #}

    #}

#}


print("OK")