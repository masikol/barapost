# -*- coding: utf-8 -*-
# This module defines functions, which are necessary for databse creating.

import os
import sys
import re
from glob import glob
from time import sleep
from threading import Thread, Event
from gzip import open as open_as_gzip
from subprocess import Popen as sp_Popen

import urllib.request
from src.lingering_https_get_request import lingering_https_get_request

from src.barapost_local_modules.barapost_spec import configure_acc_dict
from src.barapost_local_modules.related_replicons import search_for_related_replicons

import src.taxonomy as taxonomy
from src.platform import platf_depend_exit
from src.check_connection import check_connection
from src.printlog import printlog_info, printlog_info_time, printlog_error
from src.printlog import printlog_error_time, printn, getwt, log_info
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
# end if


def retrieve_fastas_by_acc(acc_dict, db_dir, local_fasta):
    # Function downloads set of records from Genbank according to accessions passed to it.
    # Downloaded FASTA file will be placed in 'db_dir' directory and named 'local_seq_set.fasta'

    # :param acc_dict: dictionary comntaining accession data of hits;
    # :type acc_dict: dict<str: tuple<str, str, int>>;
    # :param db_dir: path to directory in which downloaded FASTA file will be placed;
    # :type db_dir: str;
    # :param local_fasta: path to file with reference sequences to be included in database;
    # :type local_fasta: str;

    # Path to file with current chunk (see below "100 accession numbers...")
    tmp_fasta = os.path.join(db_dir, "tmp.fasta")

    accessions = tuple(set(acc_dict.keys()))
    if len(accessions) == 0: # just in case
        return
    # end if

    # 100 accession numbers in order not to make too long URL
    # Download genomes by chunks of 100 sequences.
    max_accnum = 100
    i = 0
    accnum = len(accessions)

    while i < accnum:

        curr_accessions = accessions[i: i + max_accnum] # slice chunk

        accs_del_comma = ','.join(curr_accessions) # accessions must be separated by comma in url
        # E-utilities provide a possibility to download records from Genbank by accessions.
        retrieve_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
db=nuccore&id={}&rettype=fasta&retmode=text".format(accs_del_comma)
        log_info("Retrieve URL: `{}`".format(retrieve_url))

        # GNU wget utility is safer, but there can be presence of absence of it :)
        wget_util = "wget"
        util_found = False
        for d in os.environ["PATH"].split(os.pathsep):
            if os.path.isdir(d) and wget_util in os.listdir(d):
                util_found = True
                break
            # end if
        # end for

        print()
        printlog_info("{} - Downloading {} reference sequences...".format(getwt(), len(curr_accessions)))

        if util_found:
            # If we have wget -- just use it

            wget_cmd = 'wget --no-check-certificate "{}" -O {}'.format(retrieve_url, tmp_fasta)
            pipe = sp_Popen(wget_cmd, shell=True)
            pipe.communicate()
            if pipe.returncode != 0:
                printlog_error_time("Error occured while downloading reference sequences")
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
                printlog_info("\r{} - {} MB downloaded ".format(getwt(), fsize))
            # end def download_waiter

            error = True
            while error:
                try:
                    waiter = Thread(target=download_waiter, args=(stop_wait,)) # create thread
                    stop_wait.set() # raise the flag
                    waiter.start() # start waiting
                    urllib.request.urlretrieve(retrieve_url, tmp_fasta) # retrieve FASTA file
                except OSError as err:
                    printlog_error_time("Error occured while downloading fasta file.")
                    printlog_error( str(err) )
                    printlog_error("`barapost-local.py` will try again in 30 seconds")
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

        printlog_info_time("Downloading is completed")

        # Write chunk to result fasta file
        with open(tmp_fasta, 'r') as infile, open(local_fasta, 'a') as outfile:
            outfile.write(infile.read())
        # end with

        # Remove temp chunk file
        os.unlink(tmp_fasta)
        i += max_accnum # go to next chunk
    # end while
# end def retrieve_fastas_by_acc


def verify_cl_accessions(accs_to_download, acc_dict):
    # Function checks existance of GenBank records that correspond to accessions
    #   specified with '-s' option. After checking the function fulills 'acc_fict'.

    # :param accs_to_download: list of accessions from command line ('-s');
    # :type accs_to_download: list<str>;
    # :param acc_dict: dictionary {<ACCESSION>: <HIT_DEFINITION>};
    # :type acc_dict: dict<str: tuple<str>>;

    check_connection("https://www.ncbi.nlm.nih.gov/")

    printlog_info_time("Verifying `-s` accessions...")
    sys.stdout.write("0/{}".format(len(accs_to_download)))

    for i, acc in enumerate(accs_to_download):

        server = "eutils.ncbi.nlm.nih.gov"
        url = "/entrez/eutils/esummary.fcgi?db=nuccore&id={}".format(acc)
        text = lingering_https_get_request(server, url, "record's name", acc)

        name = re.search(r"\<Item Name=\"Title\" Type=\"String\"\>(.+)\</Item\>", text)

        if name is None:
            printlog_info("Cannot find GenBank record with accession '{}'".format(acc))
            platf_depend_exit(1)
        else:
            name = name.group(1)
        # end if

        acc_dict[acc] = name
        sys.stdout.write("\r{}/{}".format(i+1, len(accs_to_download)))
    # end for
    print()
    printlog_info_time("OK.")
# end def verify_cl_accessions


def add_lambda_phage(local_fasta, taxonomy_path):
    # Function adds control sequence of nanopore lambda phase DNA-CS
    #    to 'local_fasta'.
    #
    # :param local_fasta: path to file with reference sequences to be included in database;
    # :type local_fasta: str;
    # :param taxonomy_path: path to taxonomy file;
    # :type taxonomy_path: str;

    print()
    printlog_info_time("Adding lambda phage control sequence...")

    # sys.path[0] is directory containing the script that was used to invoke the Python interpreter.
    # We will use it to get path to file with lambda's sequence.
    lambda_fpath = os.path.join(
        os.path.dirname(sys.path[0]),
        "lambda_control",
        "nanopore_lambda_DNA-CS_control.fasta.gz")

    # Check file existance
    if not os.path.exists(lambda_fpath):
        printlog_error_time("Error: cannot find lambda phage control sequence: '{}'".format(lambda_fpath))
        platf_depend_exit(1)
    # end if

    # Read lambda's sequence
    with open_as_gzip(lambda_fpath, 'rb') as lambda_file:
        lambda_fasta = lambda_file.read()
    # end with

    # Write it to db fasta file
    with open(local_fasta, 'wb') as db_fasta_file:
        db_fasta_file.write(lambda_fasta)
    # end with

    # Save lambda's taxonomy
    taxonomy.save_taxonomy_directly(taxonomy_path, "LAMBDA", "Lambda-phage-nanopore-control")

    printlog_info_time(" ok")
# end def add_lambda_phage


def build_local_db(tax_annot_res_dir, acc_fpath, your_own_fasta_lst, accs_to_download, use_index):
    # Function creates a database with utilities from 'blast+' toolkit
    #     according to acc_dict and your_own_fasta_lst.
    #
    # :param tax_annot_res_dir: path to current result directory
    #   (each processed file has it's own result directory);
    # :type tax_annot_res_dir: str;
    # :param acc_fpath: path to file "hits_to_download.tsv";
    # :type acc_fpath: str;
    # :param your_own_fasta_lst: list of user's fasta files to be included in database;
    # :type your_own_fasta_lst: list<str>;
    # :param accs_to_download: list of accessions from command line ('-s');
    # :type accs_to_download: list<str>;
    # :param use_index: whether to use index;
    # :type use_index: str;

    # Returns path to created database.

    # Path to directory in which database will be placed
    db_dir = os.path.join(tax_annot_res_dir, "local_database")
    # Path to DBM taxonomy file
    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy", "taxonomy.tsv")

    try:
        os.makedirs(db_dir)
    except OSError:
        #If this directory exists

        while True:
            if len(os.listdir(db_dir)) == 0:
                # If db directory is empty -- break and build a database
                break
            else:
                print()
                printlog_info("Database directory is not empty:")
                printlog_info("  `{}`".format(os.path.abspath(db_dir)))
                printlog_info("Here is it's content:")
                for i, fname in enumerate(os.listdir(os.path.abspath(db_dir))):
                    printlog_info(" {}. `{}`".format(i+1, fname))
                # end for
                reply = input("""\nPress ENTER to start classification using existing database.
Enter 'r' to remove all files in this directory and create the database from the beginning:>>""")

                if reply == "":
                    # Do not build a database, just return path to it.
                    printlog_info("You have chosen to use extant database.")

                    # Return path to DB located in this directory
                    dbpath = next(iter(os.listdir(db_dir)))
                    dbpath = dbpath.partition(".fasta")[0] + dbpath.partition(".fasta")[1] # remove all after '.fasta'

                    return os.path.join(db_dir, dbpath)

                elif reply == 'r':

                    printlog_info("You have chosen to rebuild the database.")
                    # Rename old classification files and write actual data to new one:
                    old_classif_dirs = filter(
                        lambda d: os.path.exists(os.path.join(d, "classification.tsv")),
                        glob(os.path.join(tax_annot_res_dir, "*")) )
                    old_classif_files = tuple(map(lambda f: os.path.join(f, "classification.tsv"),
                        old_classif_dirs))

                    if len(old_classif_files) > 0:
                        print()
                        printlog_info("Renaming old classification files:")
                        for classif_file in old_classif_files:
                            rename_file_verbosely(classif_file)
                        # end for
                    # end if

                    # Empty database directory
                    for file in glob("{}{}*".format(db_dir, os.sep)):
                        os.unlink(file)
                    # end for

                    # Break from the loop in order to build a database
                    break
                else:
                    print("Invalid reply: `{}`\n".format(reply))
                    continue
                # end if
            # end if
        # end while
    # end try

    # It is a dictionary of accessions and record names.
    # Accessions are keys, record names are values.
    acc_dict = configure_acc_dict(acc_fpath, your_own_fasta_lst, accs_to_download)

    if len(accs_to_download) != 0:
        verify_cl_accessions(accs_to_download, acc_dict)
    # end if

    # Retrieve already existing taxonomy data from taxonomy file
    tax_exist_accs = taxonomy.get_tax_keys(taxonomy_path)

    # If accession file does not exist and execution has reached here -- everything is OK --
    #    we are building a database from user's files only.
    if len(acc_dict) != 0:
        print()

        print("""Following sequences (and all replicons related to them)
  will be downloaded from Genbank for further taxonomic classification
  on your local machine:\n""")
        printlog_info("Following sequences (and all replicons related to them) \
will be downloaded from Genbank for further taxonomic classification \
on your local machine:")
        for i, acc in enumerate(acc_dict.keys()):
            printlog_info(" {}. {} - `{}`".format(i+1, acc, acc_dict[acc]))
        # end for

        search_for_related_replicons(acc_dict)

        printlog_info_time("Completing taxonomy file...")
        for i, acc in enumerate(acc_dict.keys()):
            if not acc in tax_exist_accs:
                taxonomy.find_taxonomy(acc, acc_dict[acc][1], taxonomy_path)
            # end if
            # Accessions can be of different length
            printn("\r{} - {}: {}/{}".format(getwt(), acc, i+1, len(acc_dict)) + " "*10 + "\b"*10)
        # end for
        print()
        printlog_info_time("Taxonomy file is consistent.")
    # end if

    local_fasta = os.path.join(db_dir, "local_seq_set.fasta") # path to downloaded FASTA file

    add_lambda_phage(local_fasta, taxonomy_path) # add lambda phage control sequence

    retrieve_fastas_by_acc(acc_dict, db_dir, local_fasta) # download main fasta data from GenBank

    # Add 'your own' fasta files to database
    if not len(your_own_fasta_lst) == 0:

        # This variable counts sequences from local files.
        # It is necessary for not allowing duplicated accessions.
        own_seq_counter = 0

        # Check if these files are assembly made by SPAdes or a5
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
                with open (local_fasta, 'a') as fasta_db:
                    for assm_path in assm_lst:
                        printlog_info("Adding `{}` to database...".format(os.path.basename(assm_path)))

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
                                    taxonomy.save_taxonomy_directly(taxonomy_path, own_acc, own_def)
                                    line = ">" + "{} {}".format(own_acc, own_def)
                                # end if
                                fasta_db.write(line + '\n')
                            # end for
                        # end with
                    # end for
                # end with
            # end if
        # end for

        with open(local_fasta, 'a') as fasta_db:
            for own_fasta_path in your_own_fasta_lst:
                printlog_info("Adding `{}` to database...".format(os.path.basename(own_fasta_path)))

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
                            taxonomy.save_taxonomy_directly(taxonomy_path, own_acc, line[1:])
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
    print()
    printn("{} - Formatting accessions...".format(getwt()))
    log_info("Formatting accessions...")
    corrected_path = os.path.join(db_dir, "corrected_seqs.fasta")
    with open(local_fasta, 'r') as source_file, open(corrected_path, 'w') as dest_file:
        for line in source_file:
            if line.startswith('>'):
                line = line.strip()
                acc, seq_name = (line.partition(' ')[0], line.partition(' ')[2])
                acc = acc.partition('.')[0]
                seq_name = seq_name.replace(' ', '_') # remove spaces
                seq_name = re.sub(r'[^\x00-\x7F]+', '_', seq_name) # remove non-ascii chars
                line = ' '.join( (acc, seq_name) ) + '\n'
            # end if
            dest_file.write(line)
        # end for
    # end with
    os.unlink(local_fasta)
    os.rename(corrected_path, local_fasta)
    sys.stdout.write("\r{} - Formatting accessions... ok".format(getwt()))
    log_info("Formatting accessions done.")

    # Configure command line
    make_db_cmd = "makeblastdb -in {} -parse_seqids -dbtype nucl".format(local_fasta)
    exit_code = os.system(make_db_cmd) # make a blast-format database
    if exit_code != 0:
        printlog_error_time("Error occured while making the database")
        platf_depend_exit(exit_code)
    # end if

    print("\033[1A{} - Database is successfully created: `{}`\n".format(getwt(), local_fasta))
    log_info("Database is successfully created: `{}`".format(local_fasta))

    if use_index == "true":
        printlog_info_time("Database index creating started")
        # Configure command line
        make_index_cmd = "makembindex -input {} -iformat blastdb -verbosity verbose".format(local_fasta)
        exit_code = os.system(make_index_cmd) # create an index for the database
        if exit_code != 0:
            printlog_info_time("Error occured while creating database index")
            platf_depend_exit(exit_code)
        # end if

        printlog_info_time("Database index has been successfully created")
    # end if

    # Gzip downloaded FASTA file
    printlog_info_time("Gzipping FASTA file: `{}`".format(local_fasta))

    if gzip_util_found:
        os.system("{} -v {}".format(gzip_util, local_fasta))
    else:
        # form .fasta.gz file 'by hand'
        with open(local_fasta, 'rb') as fasta_file, open_as_gzip(local_fasta+".gz", "wb") as fagz_file:
            shutil_copyfileobj(fasta_file, fagz_file)
        # end with
        os.unlink(local_fasta) # remove source FASTA file, not the database
    # end if

    return local_fasta
# end def build_local_db
