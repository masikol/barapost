# -*- coding: utf-8 -*-
# This module defines functions, which are necessary for databse creating.

import os
import sys
import re
import shelve
from subprocess import Popen as sp_Popen, PIPE as sp_PIPE
from time import sleep
from glob import glob
from threading import Thread, Event

import urllib.request
from src.lingering_https_get_request import lingering_https_get_request

from src.barapost_local_modules.barapost_spec import configure_acc_dict
from src.barapost_local_modules.related_replicons import search_for_related_replicons

from src.platform import platf_depend_exit
from src.check_connection import check_connection
from src.lineage import download_lineage, save_own_seq_taxonomy
from src.printlog import printl, printn, println, err_fmt, getwt
from src.filesystem import rename_file_verbosely, OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped


# GNU gzip utility is faster, but there can be presence of absence of it :)
gzip_util = "gzip"
gzip_util_found = False
for directory in os.environ["PATH"].split(os.pathsep):
    if os.path.isdir(directory) and gzip_util in os.listdir(directory):
        gzip_util_found = True
        break
    # end if
# end for

if not gzip_util_found:
    from shutil import copyfileobj as shutil_copyfileobj
    from gzip import open as open_as_gzip
# end if


def retrieve_fastas_by_acc(acc_dict, db_dir, logfile_path):
    """
    Function downloads set of records from Genbank according to accessions passed to it.
    Downloaded FASTA file will be placed in 'db_dir' directory and named 'local_seq_set.fasta'

    :param acc_dict: dictionary comntaining accession data of hits;
    :type acc_dict: dict<str: tuple<str, str, int>>;
    :param db_dir: path to directory in which downloaded FASTA file will be placed;
    :type db_dir: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;

    Returns path to downloaded FASTA file of 'str'.
    """

    local_fasta = os.path.join(db_dir, "local_seq_set.fasta") # path to downloaded FASTA file
    tmp_fasta = os.path.join(db_dir, "tmp.fasta") # path to file with current chunk (see below "100 accession numbers...")

    accessions = tuple(set(acc_dict.keys()))
    if len(accessions) == 0: # just in case
        return local_fasta
    # end if

    # 100 accession numbers in order not to make too long URL
    # Download genomes by chunks of 100 sequences.
    max_accnum = 100
    i = 0
    accnum = len(accessions)

    while i < accnum:

        curr_accessions = accessions[i: i + max_accnum] # slice chunk

        accs_del_comma = ','.join(curr_accessions) # GI numbers must be separated by comma in url
        # E-utilities provide a possibility to download records from Genbank by accessions.
        retrieve_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text".format(accs_del_comma)

        # GNU wget utility is safer, but there can be presence of absence of it :)
        wget_util = "wget"
        util_found = False
        for directory in os.environ["PATH"].split(os.pathsep):
            if os.path.isdir(directory) and wget_util in os.listdir(directory):
                util_found = True
                break
            # end if
        # end for

        printl(logfile_path, "\n{} - Downloading {} reference sequences...".format(getwt(), len(curr_accessions)))

        if util_found:
            # If we have wget -- just use it

            wget_cmd = 'wget "{}" -O {}'.format(retrieve_url, tmp_fasta)
            pipe = sp_Popen(wget_cmd, shell=True)
            pipe.communicate()
            if pipe.returncode != 0:
                printl(logfile_path, err_fmt("error while downloading reference sequences"))
                platf_depend_exit(pipe.returncode)
            # end if

        else:
            # If there are no wget -- we will download sequences with Python disposal
            
            stop_wait = Event() # a flag variable that will signal waiter-function to stop executing

            def download_waiter(stop_wait):
                """
                Function waits untill 'local_fasta' file is downloaded.
                It prints size of downloaded data to console during downloading.
                This function just waits -- it won't bring you the menu :).
                """
                # Wait untill downloading starts
                while not os.path.exists(tmp_fasta):
                    if not stop_wait.is_set():
                        return
                    # end if
                    sleep(1)
                # end while

                MB_size = 1024**2 # we will divide by it to get megabytes

                while stop_wait.is_set():
                    # Get size of downloaded data
                    fsize = round(os.path.getsize(tmp_fasta) / MB_size, 1) # get megabytes
                    printn("\r{} - {} MB downloaded ".format(getwt(), fsize))
                    sleep(1) # instant updates are not necessary
                # end while
                
                # Print total size of downloaded file (it can be deleted by this time)
                try:
                    fsize = round(os.path.getsize(tmp_fasta) / MB_size, 1)
                except OSError:
                    # We can pass this ecxeption -- we do delete this file if downloading crushes
                    # And this function just waits :)
                    pass
                # end try
                printl(logfile_path, "\r{} - {} MB downloaded ".format(getwt(), fsize))
            # end def download_waiter

            error = True
            while error:
                try:
                    waiter = Thread(target=download_waiter, args=(stop_wait,)) # create thread
                    stop_wait.set() # raise the flag
                    waiter.start() # start waiting
                    urllib.request.urlretrieve(retrieve_url, tmp_fasta) # retrieve FASTA file
                except Exception as err:
                    printl(logfile_path, err_fmt("error while downloading fasta file"))
                    printl(logfile_path,  str(err) )
                    printl(logfile_path, "'barapost-local.py' will try again in 30 seconds")
                    if os.path.exists(tmp_fasta):
                        os.unlink(tmp_fasta)
                    # end if
                    sleep(30)
                else:
                    error = False
                finally:
                    stop_wait.clear() # lower the flag
                    waiter.join() # main thread will wait until waiter function ends it's work
                # end try
            # end while
        # end if

        printl(logfile_path, "{} - Downloading is completed".format(getwt()))

        # Write chunk to result fasta file
        with open(tmp_fasta, 'r') as infile, open(local_fasta, 'a') as outfile:
            outfile.write(infile.read())
        # end with

        # Remove temp chunk file
        os.unlink(tmp_fasta)
        i += max_accnum # go to next chunk
    # end while

    return local_fasta
# end def retrieve_fastas_by_acc


def verify_cl_accessions(accs_to_download, logfile_path, acc_dict):
    """
    Function checks existance of GenBank records that correspond to accessions
      specified with '-s' option. After checking the function fulills 'acc_fict'.

    :param accs_to_download: list of accessions from command line ('-s');
    :type accs_to_download: list<str>;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    :param acc_dict: dictionary {<ACCESSION>: <HIT_DEFINITION>};
    :type acc_dict: dict<str: tuple<str>>;
    """

    check_connection("https://www.ncbi.nlm.nih.gov/")

    printl(logfile_path, "{} - Verifying '-s' accessions...".format(getwt()))
    sys.stdout.write("0/{}".format(len(accs_to_download)))

    for i, acc in enumerate(accs_to_download):

        server = "eutils.ncbi.nlm.nih.gov"
        url = "/entrez/eutils/esummary.fcgi?db=nuccore&id={}".format(acc)
        text = lingering_https_get_request(server, url, logfile_path, "record's name", acc)

        name = re.search(r"\<Item Name=\"Title\" Type=\"String\"\>(.+)\</Item\>", text)

        if name is None:
            printl(logfile_path, "Cannot find GenBank record with accession '{}'".format(acc))
            platf_depend_exit(1)
        else:
            name = name.group(1)
        # end if

        acc_dict[acc] = name
        sys.stdout.write("\r{}/{}".format(i+1, len(accs_to_download)))
    # end for
    println(logfile_path, "\n{} - OK.".format(getwt()))
# end def verify_cl_accessions


def build_local_db(tax_annot_res_dir, acc_fpath, your_own_fasta_lst, accs_to_download, logfile_path):
    """
    Function creates a database with utilities from 'blast+' toolkit
        according to acc_dict and your_own_fasta_lst.

    :param tax_annot_res_dir: path to current result directory (each processed file has it's own result directory);
    :type tax_annot_res_dir: str;
    :param acc_fpath: path to file "hits_to_download.tsv";
    :type acc_fpath: str;
    :param your_own_fasta_lst: list of user's fasta files to be included in database;
    :type your_own_fasta_lst: list<str>;
    :param accs_to_download: list of accessions from command line ('-s');
    :type accs_to_download: list<str>;
    :param logfile_path: path to logfile;
    :type logfile_path: str;

    Returns path to created database.
    """

    db_dir = os.path.join(tax_annot_res_dir, "local_database") # path to directory in which database will be placed
    # Path to DBM taxonomy file
    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy","taxonomy")

    try:
        os.makedirs(db_dir)
    except OSError as err:
        #If this directory exists

        while True:
            if len(os.listdir(db_dir)) == 0:
                # If db directory is empty -- break and build a database
                break
            else:
                printl(logfile_path, "\nDatabase exists in following directory:")
                printl(logfile_path, "  '{}'".format(os.path.abspath(db_dir)))
                reply = input("""\nPress ENTER to use existing database.
Enter 'r' to remove all files in this directory and create the database from the beginning:>>""")

                if reply == "":
                    # Do not build a database, just return path to it.
                    printl(logfile_path, "") # just print blank line

                    # Return path to DB located in this directory
                    dbpath = next(iter(os.listdir(db_dir)))
                    dbpath = dbpath.partition(".fasta")[0] + dbpath.partition(".fasta")[1] # remove all after '.fasta'

                    return os.path.join(db_dir, dbpath)
                
                elif reply == 'r':

                    # Rename old classification files and write actual data to new one:
                    old_classif_dirs = filter(
                        lambda d: os.path.exists(os.path.join(d, "classification.tsv")),
                        glob(os.path.join(tax_annot_res_dir, "*")) )
                    old_classif_files = tuple(map(lambda f: os.path.join(f, "classification.tsv"),
                        old_classif_dirs))

                    if len(old_classif_files) > 0:
                        printl(logfile_path, "\n Renaming old classification files:")
                        for classif_file in old_classif_files:
                            rename_file_verbosely(classif_file, logfile_path)
                        # end for
                    # end if

                    # Empty database directory
                    for file in glob("{}{}*".format(db_dir, os.sep)):
                        os.unlink(file)
                    # end for

                    # # Empty taxonomy file (WHY???)
                    # with shelve.open(taxonomy_path, 'n') as tax_file:
                    #     pass
                    # # end with

                    # Break from the loop in order to build a database
                    break
                else:
                    print("Invalid reply: '{}'\n".format(reply))
                    continue
                # end if
            # end if
        # end while
    # end try

    # It is a dictionary of accessions and record names.
    # Accessions are keys, record names are values.
    acc_dict = configure_acc_dict(acc_fpath, your_own_fasta_lst, logfile_path)

    if len(accs_to_download) != 0:
        verify_cl_accessions(accs_to_download, logfile_path, acc_dict)
    # end if

    # Retrieve already existing taxonomy data from taxonomy file
    with shelve.open(taxonomy_path, 'c') as tax_file:
        tax_exist_accs = tuple( tax_file.keys() )
    # end with

    # If accession file does not exist and execution has reached here -- everything is OK --
    #    we are building a database from user's files only.
    if len(acc_dict) != 0:
        print()

        printl(logfile_path, """Following sequences (and all replicons related to them)
  will be downloaded from Genbank for further taxonomic classification
  on your local machine:\n""")
        for i, acc in enumerate(acc_dict.keys()):
            printl(logfile_path, " {}. {} - '{}'".format(i+1, acc, acc_dict[acc]))
        # end for

        # Get list of GI numbers.
        search_for_related_replicons(acc_dict, logfile_path)

        printl(logfile_path, "{} - Completing taxonomy file...".format(getwt()))
        with shelve.open(taxonomy_path, 'c') as tax_file:
            for i, acc in enumerate(acc_dict.keys()):
                if not acc in tax_exist_accs:
                    download_lineage(acc, acc_dict[acc][1], tax_file, logfile_path)
                # end if
                printn("\r{} - {}: {}/{}".format(getwt(), acc, i+1, len(acc_dict)) + " "*10 + "\b"*10)
                # Accessions can be of different length
            # end for
        # end with
        printl(logfile_path, "\n{} - Taxonomy file is consistent.".format(getwt()))
    # end if

    local_fasta = retrieve_fastas_by_acc(acc_dict, db_dir, logfile_path) # download fasta file

    # Add 'your own' fasta files to database
    if not len(your_own_fasta_lst) == 0:

        # This variable counts sequences from local files.
        # It is necessary for not allowing duplicated accessions.
        own_seq_counter = 0

        # Check if these files are assembly made by SPAdese or a5
        spades_patt = r">NODE_[0-9]+" # this pattern will match sequence IDs generated y SPAdes
        spades_counter = 0 # variable counts number of SPAdes assembly files
        spades_assms = list() # this list will contain paths to SPAdes assembly files
        a5_patt = r">scaffold_[0-9]+" # this pattern will match sequence IDs generated y a5
        a5_counter = 0 # variable counts number of a5 assembly files
        a5_assms = list() # this list will contain paths to a5 assembly files

        for own_fasta_path in reversed(your_own_fasta_lst):

            how_to_open = OPEN_FUNCS[ is_gzipped(own_fasta_path) ]
            fmt_func = FORMATTING_FUNCS[ is_gzipped(own_fasta_path) ]

            with how_to_open(own_fasta_path) as fasta_file:
                first_seq_id = fmt_func(fasta_file.readline()) # get the first line in file (the first seq ID)
            # end with

            # if we've got SPAdes assembly
            if not re.search(spades_patt, first_seq_id) is None:
                spades_counter += 1
                spades_assms.append(own_fasta_path)
                continue
            # end if
            
            # if we've got a5 assembly
            if not re.search(a5_patt, first_seq_id) is None:
                a5_counter += 1
                a5_assms.append(own_fasta_path)
                continue
            # end if
        # end for

        # Include assembly files in multi-fasta file
        for counter, assm_lst in zip((spades_counter, a5_counter), (spades_assms, a5_assms)):

            # If there are more than one file with assembly of one assembler,
            #    we need to distinguish these files (i.e. these assemblies).
            # Otherwise they will be processed just like any other file
            if counter > 1:

                # Remove this files from list -- they will be processed in a specific way
                for file in assm_lst:
                    your_own_fasta_lst.remove(file)
                # end for

                # If there are any equal basenames of files, absolute paths to these files will be used in seq IDs:
                assm_basenames = list( map(os.path.basename, assm_lst) ) # get basenames

                # If this sum is equal to length of 'assm_basenames' -- there are no duplicated basenames.
                # So, there is no need to use absolute paths.
                dedpul_sum = sum( map(assm_basenames.count, assm_basenames) )

                # Path conversion according to 'deduplication sum':
                if dedpul_sum == len(assm_basenames):
                    def_convert_func = os.path.basename
                    # assm_lst = map(os.path.basename, assm_lst)
                else:
                    def_convert_func = os.path.abspath
                    # assm_lst = map(os.path.abspath, assm_lst)
                # end if

                # Add assembled sequences to database
                with open (local_fasta, 'a') as fasta_db, shelve.open(taxonomy_path, 'c') as tax_file:
                    for assm_path in assm_lst:
                        println(logfile_path, "\n{} - Adding '{}' to database...".format(getwt(), os.path.basename(assm_path)))

                        how_to_open = OPEN_FUNCS[ is_gzipped(assm_path) ]
                        fmt_func = FORMATTING_FUNCS[ is_gzipped(assm_path) ]
                        with how_to_open(assm_path) as fasta_file:
                            for line in fasta_file:
                                line = fmt_func(line)
                                # You can find comments to "OWN_SEQ..." below.
                                # Paths will be written to seq IDs in following way:
                                #   (_/some/happy/path.fastq_)
                                # in order to retrieve them securely with regex later.
                                if line.startswith('>'):
                                    own_seq_counter += 1
                                    own_acc = "OWN_SEQ_{}".format(own_seq_counter)
                                    own_def = "(_{}_)_".format(def_convert_func(assm_path)) + line[1:]
                                    own_def = own_def.replace(' ', '_')
                                    tax_file[own_acc] = own_def
                                    line = ">" + "{} {}".format(own_acc, own_def)
                                # end if
                                fasta_db.write(line + '\n')
                            # end for
                        # end with
                    # end for
                # end with
            # end if
        # end for

        with open(local_fasta, 'a') as fasta_db, shelve.open(taxonomy_path, 'c') as tax_file:
            for own_fasta_path in your_own_fasta_lst:
                println(logfile_path, "\n{} - Adding '{}' to database...".format(getwt(), os.path.basename(own_fasta_path)))

                how_to_open = OPEN_FUNCS[ is_gzipped(own_fasta_path) ]
                fmt_func = FORMATTING_FUNCS[ is_gzipped(own_fasta_path) ]
                with how_to_open(own_fasta_path) as fasta_file:
                    for line in fasta_file:
                        line = fmt_func(line)
                        # 'makeblastdb' considers first word (sep. is space) as sequence ID
                        #   and throws an error if there are duplicated IDs.
                        # In order not to allow this duplication we'll create our own sequence IDs:
                        #   'OWN_SEQ_<NUMBER>' and write it in the beginning of FASTA record name.
                        if line.startswith('>'):
                            own_seq_counter += 1
                            own_acc = "OWN_SEQ_{}".format(own_seq_counter)
                            save_own_seq_taxonomy(line[1:], own_acc, tax_file)
                            line = ">" + own_acc + ' ' + line[1:].replace(' ', '_')
                        # end if
                        fasta_db.write(line + '\n')
                    # end for
                # end with
            # end for
        # end with
    # end if

    # 'lcl|ACCESSION...' entries can be given with '.1'
    #   (or '.2', whatever) terminus by blastn.
    # There is no '.1' terminus in taxonomy file.
    # Therefore we will prune accessions in advance.
    println(logfile_path, "\n{} - Formatting accessions...".format(getwt()))
    corrected_path = os.path.join(db_dir, "corrected_seqs.fasta")
    with open(local_fasta, 'r') as source_file, open(corrected_path, 'w') as dest_file:
        for line in source_file:
            if line.startswith('>'):
                line = line.strip()
                acc, seq_name = (line.partition(' ')[0], line.partition(' ')[2])
                acc = acc.partition('.')[0]
                line = ' '.join( (acc, seq_name.replace(' ', '_')) ) + '\n'
            # end if
            dest_file.write(line)
        # end for
    # end with
    os.unlink(local_fasta)
    os.rename(corrected_path, local_fasta)
    println(logfile_path, "\r{} - Formatting accessions... ok".format(getwt()))

    # Configure command line
    make_db_cmd = "makeblastdb -in {} -parse_seqids -dbtype nucl".format(local_fasta)
    exit_code = os.system(make_db_cmd) # make a blast-format database
    if exit_code != 0:
        printl(logfile_path, err_fmt("error while making the database"))
        platf_depend_exit(exit_code)
    # end if
    
    printl(logfile_path, """\033[1A{} - Database is successfully created:
  '{}'\n""".format(getwt(), local_fasta))

    printl(logfile_path, "{} - Database index creating started".format(getwt()))
    # Configure command line
    make_index_cmd = "makembindex -input {} -iformat blastdb -verbosity verbose".format(local_fasta)
    exit_code = os.system(make_index_cmd) # create an index for the database
    if exit_code != 0:
        printl(logfile_path, err_fmt("error while creating database index"))
        platf_depend_exit(exit_code)
    # end if

    printl(logfile_path, "{} - Database index has been successfully created\n".format(getwt()))

    # Gzip downloaded FASTA file
    printl(logfile_path, "Gzipping FASTA file:\n '{}'\n".format(local_fasta))

    if gzip_util_found:
        os.system("{} {}".format(gzip_util, local_fasta))
    else:
        # form .fasta.gz file 'by hand'
        with open(local_fasta, 'rb') as fasta_file, open_as_gzip(local_fasta+".gz", "wb") as fagz_file:
            shutil_copyfileobj(fasta_file, fagz_file)
        # end with
        os.unlink(local_fasta) # remove source FASTA file, not the database
    # end if

    return local_fasta
# end def build_local_db