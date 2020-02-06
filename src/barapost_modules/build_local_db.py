# -*- coding: utf-8 -*-

from src.printlog import printl, printn, println, err_fmt, getwt
from src.platform import platf_depend_exit
from src.filesystem import rename_file_verbosely, OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped
from src.check_connection import check_connection
from src.barapost_modules.lineage import download_lineage

import urllib.request
from urllib.error import HTTPError

import os
import shelve
from glob import glob
from time import sleep
from re import search as re_search
from threading import Thread


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


def search_for_related_replicons(acc, acc_dict):
    """
    Generator finds replicons (other chromosomes or plasmids) that are related to 
      Genbank record "discodered" by prober.py.

    :param acc: accession of a record "discovered" by prober.py;
    :type acc: str;
    :param acc_dict: dictionary {<ACCESSION>: (<GI_NUMBER>, <HIT_DEFINITION>)};
    :type acc_dict: dict<str: tuple<str>>;

    Yields tuples of a following structure:
        (<ACCESSION>, <GI_NUMBER>, <RECORD_DEFINITION>)
    """

    # Get the smallest web page that contains BioSample:
    summary_url = "https://www.ncbi.nlm.nih.gov/nuccore/{}?report=docsum&log$=seqview".format(acc)
    summary_html = urllib.request.urlopen(summary_url).read().decode("utf-8")
    # Get reference to BioSample web page:
    biosample_regex = r"href=\"/(biosample\?LinkName=nuccore_biosample&amp;from_uid=[0-9]+)"
    biosample_ref = "https://www.ncbi.nlm.nih.gov/" + re_search(biosample_regex, summary_html).group(1)
    del summary_html # let it go

    # Get BioSample web page:
    biosample_html = urllib.request.urlopen(biosample_ref).read().decode("utf-8")
    # Get reference to list nucleotide links:
    nucl_regex = r"href=\"/(nuccore\?LinkName=biosample_nuccore&amp;from_uid=[0-9]+)"
    nucl_ref = "https://www.ncbi.nlm.nih.gov/" + re_search(nucl_regex, biosample_html).group(1)
    del biosample_html # let it go

    # Get list nucleotide links:
    nucl_html = urllib.request.urlopen(nucl_ref).read().decode("utf-8")

    num_links = 0 # number of nucleotide links of on this page

    # Count these links:
    while 'ordinalpos={}">'.format(num_links+1) in nucl_html:
        num_links += 1
    # end while

    # Duplicated records (e.g. from RefSeq) should not be allowed.
    # This is the list of already encountered record definitions]
    #   (e.g. "Erwinia amylovora strain E-2 plasmid pEa-E-2, complete sequence"):
    def_list = [ acc_dict[acc][1] ] # "discovered" one is already here

    for ordinalpos in range(1, num_links+1):

        # Get text containing record definition, accession and GI number:
        payload_regex = r"(ordinalpos={}\"\>.*?GI: \</dt\>\<dd\>[0-9]+)".format(ordinalpos)
        text = re_search(payload_regex, nucl_html).group(1)

        # Get definition of a record
        definition = re_search(r"ordinalpos={}\"\>(.+)\</a\>".format(ordinalpos), text).group(1)

        if not definition in def_list:

            def_list.append(definition)

            # Get accession and GI number:
            rel_acc = re_search(r"Accession: \</dt\>\<dd>(.*?)\.", text).group(1)
            rel_gi = re_search(r"GI: \</dt\>\<dd\>([0-9]+)", text).group(1)

            acc_dict[rel_acc] = (rel_gi, definition) # update acc_dict

            yield rel_acc, rel_gi, definition
        # end if
    # end for

# end def search_for_related_replicons


def get_gi_by_acc(acc_dict, logfile_path):
    """
    Function returns GI number that corresponds to accessing passed to it.
    """

    printl(logfile_path, "\n{} - Searching for related replicons...".format(getwt()))

    gi_list = list()

    start_accs = tuple(acc_dict.keys())

    for i, acc in enumerate(start_accs):
        if not acc_dict[acc][0] in gi_list:
            try:
                gi_list.append(acc_dict[acc][0])
            except KeyError:
                printl(logfile_path, err_fmt("GI number error. Please, contact the developer"))
                platf_depend_exit(1)
            # end try
        # end if

        # Search for related replicons:
        try:
            related_repls = search_for_related_replicons(acc, acc_dict)
        except AttributeError:
            print("\nParsing error: cannot find replicons related to {}.".format(acc))
            print("Please, contact the developer")
        else:
            for rel_acc, rel_gi, definition in related_repls:
                if not rel_gi in gi_list:
                    gi_list.append(rel_gi)
                    printl(logfile_path, "\r {} - {}".format(rel_acc, definition))
                    printn("{} - Searching for related replicons...".format(getwt()))
                    # print(gi_list)
                # end if
            # end for
        # end try
    # end for

    if len(start_accs) != len(acc_dict):
        printl(logfile_path, "\r{} - {} related replicons have been found.".format(getwt(),
            len(acc_dict) - len(start_accs)))
    else:
        printl(logfile_path, "\r{} - No related replicons found.".format(getwt()))
    # end if

    return gi_list
# end def get_gi_by_acc


def retrieve_fastas_by_gi(gi_list, db_dir, logfile_path):
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
    if len(gi_list) == 0:
        return local_fasta
    # end if

    gis_del_comma = ','.join(gi_list) # GI numbers must be separated by comma in url
    # E-utilities provide us with possibility of downloading records from Genbank by GI numbers.
    retrieve_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text".format(gis_del_comma)
    
    global stop_wait # a flag variable that will signal waiter-function to stop executing

    def download_waiter():
        """
        Function waits untill 'local_fasta' file is downloaded.
        It prints size of downloaded data to console during downloading.
        This function just waits -- it won't bring you the menu :).
        """
        # Wait untill downloading starts
        while not os.path.exists(local_fasta):
            if stop_wait:
                return
            # end if
            sleep(1)
        # end while

        while not stop_wait:
            # Get size of downloaded data
            fsize = round(os.path.getsize(local_fasta) / (1024**2), 1) # get megabytes
            printn("\r{} - {} MB downloaded ".format(getwt(), fsize))
            sleep(1) # instant updates are not necessary
        # end while
        
        # Print total size of downloaded file (it can be deleted by this time)
        try:
            fsize = round(os.path.getsize(local_fasta) / (1024**2), 1)
        except OSError:
            pass # we can pass this ecxeption -- we do delete this file if downloading crushes
        # end try
        printl(logfile_path, "\r{} - {} MB downloaded \n".format(getwt(), fsize))
    # end def download_waiter

    error = True
    while error:
        try:
            waiter = Thread(target=download_waiter) # create thread
            stop_wait = False # raise the flag
            waiter.start() # start waiting
            printl(logfile_path, "\n{} - Downloading sequences for local database building started".format(getwt()))
            urllib.request.urlretrieve(retrieve_url, local_fasta) # retrieve FASTA file
        except Exception as err:
            stop_wait = True
            printl(logfile_path, err_fmt("error while downloading FASTA files"))
            printl(logfile_path,  str(err) )
            printl(logfile_path, "'barapost.py' will try again in 30 seconds")
            if os.path.exists(local_fasta):
                os.unlink(local_fasta)
            # end if
            sleep(30)
        else:
            error = False
        finally:
            stop_wait = True # lower the flag
            waiter.join() # main thread will wait until waiter function ends it's work
        # end try
    # end while

    println(logfile_path, "{} - Downloading is completed".format(getwt()))

    return local_fasta
# end def retrieve_fastas_by_gi


def build_local_db(acc_dict, tax_annot_res_dir, acc_fpath, your_own_fasta_lst, logfile_path):
    """
    Function builds a local indexed database with utilities from 'blast+' toolkit.

    :param acc_dict: a dictionary of accessions and record names
        Accession are keys, record names are values;
    :type acc_dict: dict<str, str>;
    :param tax_annot_res_dir: path to current result directory (each processed file has it's own result directory);
    :type tax_annot_res_dir: str;

    Returns path to builded local indexed database.
    """

    db_dir = os.path.join(tax_annot_res_dir, "local_database") # path to directory in which database will be placed
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
                reply = input("""\nPress ENTER to continue aligning using this database.
Enter 'r' to remove all files in this directory and build the database from the beginning:>>""")

                if reply == "":
                    # Do not build a database, just return path to it.
                    printl(logfile_path, ) # just print blank line
                    return os.path.join(db_dir, "local_seq_set.fasta")
                
                elif reply == 'r':
                    if acc_fpath is None and len(your_own_fasta_lst) == 0:
                        printl(logfile_path, err_fmt("missing data to build a database from"))
                        printl(logfile_path, """ There is no accession file in directory '{}'
 and no FASTA file have been specified with '-l' option.""")
                        platf_depend_exit(1)
                    # end if

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

                    # Empty this directory and break from the loop in order to build a database.
                    for file in glob("{}{}*".format(db_dir, os.sep)):
                        os.unlink(file)
                    # end for

                    break
                else:
                    # Ask again
                    continue
                # end if
            # end if
        # end while
    # end try

    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy","taxonomy")

    with shelve.open(taxonomy_path, 'c') as tax_file:
        tax_exist_accs = tuple( tax_file.keys() )
    # end with

    # If accession file does not exist and execution has reached here -- everything is OK --
    #    we are building local database from local files only.
    if not acc_fpath is None:
        print()
        check_connection("https://www.ncbi.nlm.nih.gov/")

        printl(logfile_path, """\nFollowing sequences (and all replicons related to them)
  will be downloaded from Genbank for further taxonomic classification
  on your local machine:\n""")
        for i, acc in enumerate(acc_dict.keys()):
            printl(logfile_path, " {}. {} - '{}'".format(i+1, acc, acc_dict[acc][1]))
        # end for

        # Get list of GI numbers.
        gi_list = get_gi_by_acc(acc_dict, logfile_path)

        printl(logfile_path, "\n{} - Completing taxonomy file...".format(getwt()))
        for acc in acc_dict.keys():
            if not acc in tax_exist_accs:
                download_lineage(acc_dict[acc][0], acc_dict[acc][1], acc, tax_annot_res_dir)
            # end if
        # end for
        printl(logfile_path, "{} - Taxonomy file is consistent.".format(getwt()))
    # end if

    local_fasta = retrieve_fastas_by_gi(gi_list, db_dir, logfile_path) # download FASTA file
    print('\n')

    # Add 'your own' FASTA files to database
    if not len(your_own_fasta_lst) == 0:

        # This variable counts sequences from local files.
        # It is necessary for accession deduplication.
        own_seq_counter = 0

        # Check if these files are SPAdes of a5 assembly
        spades_patt = r"NODE_[0-9]+" # this pattern will match sequence IDs generated y SPAdes
        spades_counter = 0 # variable counts number of SPAdes assembly files
        spades_assms = list() # this list will contain paths to SPAdes assembly files
        a5_patt = r"scaffold_[0-9]+" # this pattern will match sequence IDs generated y a5
        a5_counter = 0 # variable counts number of a5 assembly files
        a5_assms = list() # this list will contain paths to a5 assembly files

        for own_fasta_path in your_own_fasta_lst:

            how_to_open = OPEN_FUNCS[ is_gzipped(own_fasta_path) ]
            fmt_func = FORMATTING_FUNCS[ is_gzipped(own_fasta_path) ]

            with how_to_open(own_fasta_path) as fasta_file:
                first_seq_id = fmt_func(fasta_file.readline()) # get the first line in file (the first seq ID)
            # end with

            # if we've got SPAdes assembly
            if not re_search(spades_patt, first_seq_id) is None:
                spades_counter += 1
                spades_assms.append(own_fasta_path)
                continue
            # end if
            
            # if we've got SPAdes assembly
            if not re_search(a5_patt, first_seq_id) is None:
                a5_counter += 1
                a5_assms.append(own_fasta_path)
                continue
            # end if
        # end for

        for counter, assm_lst in zip((spades_counter, a5_counter), (spades_assms, a5_assms)):

            # If there are more than one file with assembly of one assembler,
            #    we need to distinguish these files (i.e. these assemblies).
            if counter > 1:

                # Remove this files from list -- they will be processed in a specific way
                for file in assm_lst:
                    your_own_fasta_lst.remove(file)
                # end for

                # If there are any equal basenames of files, absolute paths to these files will be used in seq IDs:
                assm_basenames = list( map(os.path.basename, assm_lst) ) # get basenames

                # If this sum is equal to length of 'assm_basenames' -- there are no duplicated basenames.
                # So, there is no need to use absolute paths.# Absolute path will be used otherwise.
                dedpul_sum = sum( map(assm_basenames.count, assm_basenames) )

                # Path conversion according to 'deduplication sum':
                if dedpul_sum == len(assm_basenames):
                    assm_lst = list( map(os.path.basename, assm_lst) )
                else:
                    assm_lst = list( map(os.path.abspath, assm_lst) )
                # end if

                # Add assembled sequences to database
                fasta_db = open(local_fasta, 'a')
                for assm_path in assm_lst:
                    printl(logfile_path, "{} - Adding '{}' to database...".format(getwt(), os.path.basename(assm_path)))

                    how_to_open = OPEN_FUNCS[ is_gzipped(assm_path) ]
                    with how_to_open(assm_path) as fasta_file:
                        for line in fasta_file:
                            line = line.strip()
                            # You can find comments to "OWN_SEQ..." below. I don't want to duplicate them.
                            # Paths will be written to seq IDs in following way:
                            #   (_/some/happy/path.fastq_)
                            # in order to retrieve them securely with regex later.
                            if line.startswith('>'):
                                own_seq_counter += 1
                                own_acc = "OWN_SEQ_{}".format(own_seq_counter)
                                own_def = "(_{}_)_" + line[1:]
                                with shelve.open(taxonomy_path, 'c') as tax_file:
                                    tax_file[own_acc] = own_def
                                # end with
                                line = ">" + "{} {}".format(own_acc, own_def)
                            # end if
                            fasta_db.write(line + '\n')
                        # end for
                    # end with
                # end for
                fasta_db.close()
            # end if
        # end for

        # No 'with open' here in order not to indent too much.
        fasta_db = open(local_fasta, 'a')
        for own_fasta_path in your_own_fasta_lst:
            printl(logfile_path, "{} - Adding '{}' to database...".format(getwt(), os.path.basename(own_fasta_path)))

            how_to_open = OPEN_FUNCS[ is_gzipped(own_fasta_path) ]
            with how_to_open(own_fasta_path) as fasta_file:
                for line in fasta_file:
                    line = fmt_func(line)
                    # 'makeblastdb' considers first word (sep. is space) as sequence ID
                    #   and throws an error if there are duplicate IDs.
                    # In order not to allow this duplication we'll create our own sequence IDs:
                    #   'OWN_SEQ_<NUMBER>' and write it in the beginning of FASTA record name.
                    if line.startswith('>'):
                        own_seq_counter += 1
                        own_acc = "OWN_SEQ_{}".format(own_seq_counter)
                        with shelve.open(taxonomy_path, 'c') as tax_file:
                            tax_file[own_acc] = line[1:]
                        # end with
                        line = ">" + own_acc + ' ' + line[1:]
                    # end if
                    
                    fasta_db.write(line + '\n')
                # end for
            # end with
        # end for
        fasta_db.close()
    # end if

    # Configure command line
    make_db_cmd = "makeblastdb -in {} -parse_seqids -dbtype nucl".format(local_fasta)
    exit_code = os.system(make_db_cmd) # make a blast-format database
    if exit_code != 0:
        printl(logfile_path, err_fmt("error while making the database"))
        platf_depend_exit(exit_code)
    # end if
    
    printl(logfile_path, """{} - Database is successfully created:
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