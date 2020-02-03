#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.15.b"
# Year, month, day
__last_update_date__ = "2020-01-29"

# |===== Check python interpreter version =====|

import sys

if sys.version_info.major < 3:
    print( "\nYour python interpreter version is " + "%d.%d" % (sys.version_info.major,
        sys.version_info.minor) )
    print("   Please, use Python 3.\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith("win"):
        raw_input("Press ENTER to exit:")
    # end if
    exit(1)
# end if

from src.platform import platf_depend_exit

# First search for information-providing options:

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:

    print("\n  prober.py\n  Version {}; {} edition;\n".format(__version__, __last_update_date__))
    print("DESCRIPTION:\n")
    print("""This script is designed for determinating the taxonomic position
  of nucleotide sequences by sending each of them to NCBI BLAST server and regarding the best hit.\n""")

    if "--help" in sys.argv[1:]:
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

    if "--help" in sys.argv[1:]:
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

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

import os
from re import search as re_search
from glob import glob
from src.printlog import getwt, get_full_time, printl, err_fmt

import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvd:o:p:a:g:b:e:x:",
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
    if not is_fq_or_fa(arg):
        print(err_fmt("invalid positional argument: '{}'".format(arg)))
        print("Only FAST(A/Q) files can be specified without an option in command line.")
        platf_depend_exit(1)
    # end if
    if not os.path.exists(arg):
        print(err_fmt("File '{}' does not exist!".format(arg)))
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
        if arg != "" and re_search(r"^.+@.+\..+$", arg) is None:
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
            print(err_fmt("no FASTQ or FASTA files to process found."))
            platf_depend_exit(1)
        # end if
    # end if
# end if

# Sort list of files that will be processed -- process them in alphabetical order.
fq_fa_list.sort()

from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS
from src.filesystem import remove_tmp_files, create_result_directory
from src.filesystem import is_gzipped, is_fastq, is_fasta
from src.platform import get_logfile_path

# |=== Check if there are enough sequeneces in files (>= probing_batch_size) ===|
seqs_at_all = 0
print() # just print new line
for file in fq_fa_list:
    how_to_open = OPEN_FUNCS[ is_gzipped(file) ]
    if is_fastq(file):
        seqs_at_all += sum(1 for line in how_to_open(file)) // 4  # 4 lines per record
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
if seqs_at_all < probing_batch_size:
    if ("-b" in sys.argv[1:] or "--probing_batch_size" in sys.argv[1:]):
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
    else:
        probing_batch_size = seqs_at_all
    # end if
# end if

if len(tuple(filter(is_fastq, fq_fa_list))) != 0:
    from src.fastq import fastq_packets
# end if

if len(tuple(filter(is_fasta, fq_fa_list))) != 0:
    from src.fasta import fasta_packets
# end if

packet_size = min(packet_size, probing_batch_size)

if not os.path.isdir(outdir_path):
    try:
        os.makedirs(outdir_path)
    except OSError as oserr:
        print(err_fmt("unable to create result directory '{}'".format(outdir)))
        print( str(oserr) )
        platf_depend_exit(1)
    # end try
# end if

acc_fpath = os.path.join(outdir_path, "hits_to_download.tsv")

taxonomy_dir = os.path.join(outdir_path, "taxonomy")
if not os.path.isdir(taxonomy_dir):
    os.makedirs(taxonomy_dir)
# end if
taxonomy_path = os.path.join(taxonomy_dir, "taxonomy")

from src.check_connection import check_connection

check_connection("https://blast.ncbi.nlm.nih.gov")

logfile_path = get_logfile_path("prober", outdir_path)

printl(logfile_path, "\n |=== prober.py (version {}) ===|\n".format(__version__))
printl(logfile_path, get_full_time() + "- Start working\n")

from src.prober_modules.networking import verify_taxids

organisms = verify_taxids(taxid_list, logfile_path)

from src.prober_modules.networking import lingering_https_get_request
from src.prober_modules.networking import configure_request, send_request, wait_for_align
from src.write_classification import write_classification
from src.prober_modules.prober_spec import ask_for_resumption, look_around
from src.prober_modules.prober_spec import parse_align_results_xml, write_hits_to_download
from src.prober_modules.lineage import get_lineage


#                       |===== Proceed =====|

# 1. 'previous_data' is a dict of the following structure:
#    {
#        "RID": saved_RID (str),
#        "tsv_respath": path_to_tsv_file_from_previous_run (str),
#        "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file (int),
#        "tmp_fpath": path_to_pemporary_file (str),
#        "decr_pb": value to subtract from probing_batch_size
#                  in orde to resume correctly (int)
#    }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. 'packet' is a dict of the following structure:
#    {
#        "fasta": FASTA_data_containing_query_sequences (str),
#        "qual": dictionary {seq_id: quality}
#    }
#    quality is '-' for fasta files

#                   |===== Kernel loop =====|


printl(logfile_path, " - Output directory: '{}';".format(outdir_path))
printl(logfile_path, " - Logging to '{}'".format(logfile_path))
if user_email != "":
    printl(logfile_path, " - Your email: {}".format(user_email))
# end if
printl(logfile_path, " - Probing batch size: {} sequences;".format(probing_batch_size))
printl(logfile_path, " - Packet size: {} sequences;".format(packet_size))
if not max_seq_len is None:
    printl(logfile_path, " - Maximum length of a sequence sent: {} bp;".format(max_seq_len))
# end if
printl(logfile_path, " - BLAST algorithm: {};".format(blast_algorithm))
printl(logfile_path, " - Database: nt;")
if len(organisms) > 0:
    for db_slice in organisms:
        printl(logfile_path, "   {};".format(db_slice))
    # end for
# end if

s_letter = '' if len(fq_fa_list) == 1 else 's'
printl(logfile_path, "\n {} file{} will be processed.".format( len(fq_fa_list), s_letter))
with open(logfile_path, 'a') as logfile:
    logfile.write("Here they are:\n")
    for i, path in enumerate(fq_fa_list):
        logfile.write("    {}. '{}'\n".format(i+1, path))
    # end for
# end with

printl(logfile_path, '-'*30)


# Dictionary of accessions and record names.
# Accessions are keys, record names are values.
# This dictionary is filled while processing and at the beginning of resumption.
acc_dict = dict()

# Counter of sequences processed during current run
seqs_processed = 0
# Counter of sequences processed concerning putative previous run(s)
glob_seqs_processed = 0

# Varible for stopping execution when probing batch is processed completely.
stop = False

# Iterate through found source FASTQ and FASTA files
for i, fq_fa_path in enumerate(fq_fa_list):

    # Create the result directory with the name of FASTQ of FASTA file being processed:
    new_dpath = create_result_directory(fq_fa_path, outdir_path)

    # "hname" means human readable name (i.e. without file path and extention)
    infile_hname = os.path.basename(fq_fa_path)
    infile_hname = re_search(r"(.+)\.(m)?f(ast)?(a|q)(\.gz)?$", infile_hname).group(1)

    # Look around and ckeck if there are results of previous runs of this script
    # If 'look_around' is None -- there is no data from previous run
    previous_data = look_around(outdir_path, new_dpath, fq_fa_path,
        blast_algorithm, acc_dict, probing_batch_size, logfile_path)

    if previous_data is None: # If there is no data from previous run
        num_done_seqs = 0 # number of successfully processed sequences
        tsv_res_path = "{}.tsv".format(os.path.join(new_dpath,
            "classification")) # form result tsv file path
        tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath,
            infile_hname), blast_algorithm) # form temporary file path
        saved_RID = None
    else: # if there is data from previous run
        num_done_seqs = previous_data["n_done_reads"] # get number of successfully processed sequences
        tsv_res_path = previous_data["tsv_respath"] # result tsv file sholud be the same as during previous run
        tmp_fpath = previous_data["tmp_fpath"] # temporary file sholud be the same as during previous run
        saved_RID = previous_data["RID"] # having this RID we can try to get response for last request
        probing_batch_size -= previous_data["decr_pb"]
    # end if

    glob_seqs_processed += num_done_seqs

    # Calculate total number of packets meant to be sent from current FASTA file
    packs_at_all = probing_batch_size // packet_size

    if probing_batch_size % packet_size > 0: # And this is ceiling
        packs_at_all += 1
    # end if

    pack_to_send = 1 # number of packet meant to be sent now

    packet_generator = fastq_packets if is_fastq(fq_fa_path) else fasta_packets

    # Iterate over packets left to send
    for packet in packet_generator(fq_fa_path, packet_size, num_done_seqs, max_seq_len):

        if packet["fasta"] == "":
            break
        # end if

        # Assumption that we need to send current packet
        send = True
        align_xml_text = None

        # If current packet has been already send, we can try to get it and not to send again
        if not saved_RID is None:

            resume_rtoe = 0 # we will not sleep at the very beginning of resumption

            align_xml_text = wait_for_align(saved_RID, resume_rtoe,
                pack_to_send, packs_at_all, os.path.basename(fq_fa_path), logfile_path) # get BLAST XML response
            saved_RID = None

            if align_xml_text is None:
                # Resend the packet ('send' remains True)
                pass
            elif align_xml_text == "[blastsrv4.REAL]":
                # Halve sequences and resend the packet ('send' remains True)
                packet["fasta"] = prune_seqs(packet["fasta"], 'f', 0.5)
                align_xml_text = None
            else:
                # OK -- omit 'if send' statement and write results
                send = False
            # end if
        # end if

        if send:

            printl(logfile_path, "\nGoing to BLAST (" + blast_algorithm + ")")
            printl(logfile_path, "Request number {} out of {}. Sending {} sequences.".format(pack_to_send,
                packs_at_all, len(packet["qual"])))

            while align_xml_text is None: # until successfull attempt

                request = configure_request(packet["fasta"], blast_algorithm, organisms, user_email) # get the request

                # Send the request and get BLAST XML response
                # 'align_xml_text' will be None if NCBI BLAST server rejects the request due to too large amount of data in it.

                align_xml_text = send_request(request, pack_to_send, packs_at_all,
                    os.path.basename(fq_fa_path), tmp_fpath, logfile_path)

                # If NCBI BLAST server rejects the request due to too large amount of data in it --
                #    shorten all sequences in packet twofold and resend it.
                if align_xml_text == "[blastsrv4.REAL]":
                    packet["fasta"] = prune_seqs(packet["fasta"], 'f', 0.5)
                    align_xml_text = None
                # end if
            # end while
        # end if

        # Get result tsv lines
        result_tsv_lines = parse_align_results_xml(align_xml_text,
            packet["qual"], acc_dict, logfile_path, taxonomy_path) # get result tsv lines

        # Write classification to TSV file
        write_classification(result_tsv_lines, tsv_res_path)
        # Write accessions and GI numbers of hits to TSV file
        write_hits_to_download(acc_dict, acc_fpath)

        seqs_processed += len( packet["qual"] )
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


from src.get_undr_sep_number import get_undr_sep_number

glob_seqs_processed += seqs_processed
str_about_prev_runs = ", including previous run(s)" if glob_seqs_processed > seqs_processed else ""

printl(logfile_path, '-'*20+'\n')
printl(logfile_path, " {} sequences have been processed{}\n".format(get_undr_sep_number(glob_seqs_processed),
    str_about_prev_runs))

printl(logfile_path, "Here are Genbank records that can be used for further sorting by 'barapost.py'.")
printl(logfile_path, "They are sorted by their occurence in probing batch:")

# Print accessions and record names sorted by occurence
# "-x[1][2]:": minus because we need descending order, [1] -- get tuple of "other information",
#   [2] -- get 3-rd element (occurence)
for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][2]):
    s_letter = "s" if other_info[2] > 1 else ""
    printl(logfile_path, " {} hit{} - {}, '{}'".format(get_undr_sep_number(other_info[2]),
        s_letter, acc, other_info[1]))
# end for

# Print number of unkmown sequences, if there are any:
unkn_num = glob_seqs_processed - sum( map(lambda x: x[2], acc_dict.values()) )
if unkn_num > 0:
    s_letter = "s" if unkn_num > 1 else ""
    printl(logfile_path, " {} sequence{} - No significant similarity found".format(unkn_num, s_letter))
# end if

printl(logfile_path, """\nThey are saved in following file:
    '{}'""".format(acc_fpath))
printl(logfile_path, """\nYou can edit this file before running 'barapost.py' in order to
  modify list of sequences that will be downloaded from Genbank
  and used as local (i.e. on your local computer) database by 'barapost.py'.""")

printl(logfile_path, '\n' + get_full_time() + "- Task is completed\n")
platf_depend_exit(0)