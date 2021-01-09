#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.22.a"
# Year, month, day
__last_update_date__ = "2021-01-09"

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
    sys.exit(1)
else:
    print("\nPython {}.{}.{}".format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
# end if

from src.platform import platf_depend_exit

# Firstly search for information-providing options (-h and -v):

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:

    print("\n  barapost-prober.py\n  Version {}; {} edition;\n".format(__version__, __last_update_date__))
    print("DESCRIPTION:\n")
    print("""This script is designed for taxonomic classification
  of nucleotide sequences by finding the most similar sequence
  in NCBI Nucleotide database using BLAST web service.\n""")

    if "--help" in sys.argv[1:]: # print more detailed help message
        print("""Processing entire data set with "barapost-prober.py" is not expedient,
  since servers are in public use and handling a request can linger for a long time.""")
        print("""Therefore the main goal of this script is to submit a relatively small
  probing batch (see `-b` option) of sequences to NCBI BLAST service and discover,
  what Genbank records can be downloaded and used as reference sequences in a database
  stored on local machine. Further classification will be performed by "barapost-local.py" on local machine.""")
        print("This script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.")
        print("----------------------------------------------------------\n")
        print("Default parameters:\n")
        print("""- if no input files are specified, "barapost-prober.py" processes
   all FASTQ and FASTA files in working directory;""")
        print("""- algorithm (see '-a' option): '0' (zero, i.e. megaBlast);""")
        print("""- packet forming mode (see '-c' option): '0' (zero);""")
        print("""- packet size: (see '-p' option):
   100 sequences for packet forming mode '0',
   20,000 base pairs for packet forming mode '1';""")
        print("""- probing batch size (see '-b' option): 200 sequences;""")
        print("""- database slices (see '-g' option): whole 'nr/nt' database, i.e. no slices;""")
        print("""- output directory ('-o' option): directory named 'barapost_result'
   nested in working directory;""")
        print("""- prober sends no email information ('-e' option) to NCBI;""")
        print("""- prober sends sequences intact
   (i.e. does not prune them before submission (see '-x' option));""")
        print("----------------------------------------------------------\n")
    # end if

    print("""Files that you want 'barapost-prober.py' to process should be specified as
  positional arguments (see EXAMPLE #2 in detailed (--help) help message).\n""")
    print("OPTIONS:\n")
    print("""-h (--help) --- show help message.
   '-h' -- brief, '--help' -- full;\n""")
    print("-v (--version) --- show version;\n")
    print("""-d (--indir) --- directory which contains FASTQ of FASTA files meant to be processed.
   I.e. all FASTQ and FASTA files in this direcory will be processed;\n""")
    print("""-o (--outdir) --- output directory.
   Default value: 'barapost_result';\n""")
    print("""-c (--packet-mode) --- indicates how prober calculates size of packet.
   Values:
   0 -- packet size equals number of sequences in it [default mode];
   1 -- packet size equals sum of lengths of sequences in it;
   Mode 1 performs better if input sequences are sorted by length within input file;\n""")
    print("""-p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
   It means that "barapost-prober.py" can preform multiple requests during single run.
   Value: positive integer number.
   Default value:
   100 (sequences) for mode 0;
   20,000 (base pairs) for mode 1;\n""")
    print("""-a (--algorithm) --- BLASTn algorithm to use for alignment.
   Available values: 0 for megaBlast, 1 for discoMegablast, 2 for blastn.
   Default is 0 (megaBlast);\n""")
    print("""-g (--organisms) --- TaxIDs of organisms to align your sequences against. I.e. 'nr/nt' database slices.
   Functionality of this option is totally equal to "Organism" text boxes on this BLASTn page:
    'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'.
   Format of value (see EXAMPLES #4 and #5 below.):
     <organism1_taxid>,<organism2_taxid>...
   Spaces are NOT allowed.
   Default value is full 'nr/nt' database, i.e. no slices.
   You can find your Taxonomy IDs here: 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi';\n""")
    print("""-b (--probing-batch-size) --- total number of sequences that will be submitted to BLAST server
   during 'barapost-prober.py' run.
   You can specify '-b all' to process all your sequeces by 'barapost-prober.py'.
   Value: positive integer number.
   Default value is 200;\n""")
    print("""-e (--email) --- your email. Please, specify your email when you run "barapost-prober.py",
   so that the NCBI can contact you if there is a problem. See Example #2 below;\n""")
    print("""-x (--max-seq-len) --- maximum length of a sequence that prober subits to NCBI BLAST service.
   It means that prober can prune your sequences before submission in order to spare NCBI servers.
   This feature is disabled by default;""")

    if "--help" in sys.argv[1:]:
        print("----------------------------------------------------------\n")
        print("EXAMPLES:\n")
        print("""1. Process all FASTA and FASTQ files in working directory with default settings:\n
   barapost-prober.py\n""")
        print("""2. Process all files starting with "some_fasta" in the working directory with default settings.
   Provide NCBI with your email:\n
   barapost-prober.py some_fasta* -e my.email@sth.com\n""")
        print("""3. Process one file with default settings:\n
   barapost-prober.py reads.fastq\n""")
        print("""4. Process a FASTQ file and a FASTA file with discoMegablast, packet size of 100 sequences.
   Search only among Erwinia sequences (551 is Erwinia taxid):\n
   barapost-prober.py reads_1.fastq.gz some_sequences.fasta -a discoMegablast -p 100 -g 551\n""")
        print("""5. Process all FASTQ and FASTA files in directory named `some_dir`.
   Process 300 sequences, packet size is 100 sequnces (3 packets will be sent).
   Search only among Escherichia (taxid 561) and viral (taxid 10239) sequences:\n
   barapost-prober.py -d some_dir -g 561,10239 -o outdir -b 300 -p 100\n""")
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
from src.printlog import get_full_time, printl, err_fmt

# Get command line arguments
import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvd:o:p:c:a:g:b:e:x:",
        ["help", "version", "indir=", "outdir=", "packet-size=", "packet-mode=",
        "algorithm=", "organisms=", "probing-batch-size=", "email=",
        "max-seq-len="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

is_fq_or_fa = lambda f: True if not re_search(r".*\.(m)?f(ast)?(a|q)(\.gz)?$", f) is None else False

# Default values:
fq_fa_list = list() # list of paths to file meant to be processed
indir_path = None # path to '-d' directory
outdir_path = "barapost_result"
packet_size = 100
probing_batch_size = 200
send_all = False # it will be True if '-b all' is specified
blast_algorithm = "megaBlast"
taxid_list = list() # list of TaxIDs to perform database slices
user_email = ""
max_seq_len = float("inf") # maximum length of a sequence sent to NCBI
packet_mode = 0 # mode of packet forming. 'numseqs' is default

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

# Handle command line options and values passed with them:
for opt, arg in opts:

    if opt in ("-o", "--outdir"):
        outdir_path = os.path.abspath(arg)

    elif opt in ("-p", "--packet-size"):
        try:
            packet_size = int(arg)
            if packet_size <= 0:
                raise ValueError
            # end if
        except ValueError:
            print(err_fmt("packet size (-p option) must be integer number > 1"))
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-c", "--packet-mode"):
        if arg == "0":
            packet_mode = 0
        elif arg == "1":
            packet_mode = 1
            if not "-p" in sys.argv[1:] and not "--packet-size" in sys.argv[1:]:
                packet_size = 20000 # default for mode 1
            # end if
        else:
            print(err_fmt("packet forming mode (-c option) must be 0 or 1."))
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end if

    elif opt in ("-g", "--organisms"):

        taxid_list = arg.strip().split(',')

        try: # primarily verify TaxIds
            for taxid in taxid_list:
                buff_var = int(taxid)
                if buff_var < 0:
                    raise ValueError
                # end if
            # end for
        except ValueError:
            print(err_fmt("TaxID should be a positive integer number!\a"))
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-b", "--probing-batch-size"):
        # Switch 'send_all' to True in order to process all sequences
        if arg == "all":
            send_all = True
            probing_batch_size = float("inf")
            continue
        # end if
        try:
            probing_batch_size = int(arg)
            if probing_batch_size <= 0:
                raise ValueError
            # end if
        except ValueError:
            print(err_fmt("probing batch size ('-b' option) must be positive integer number!"))
            print("Your value: {}".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-e", "--email"):
        if re_search(r"^.+@.+\..+$", arg) is None:
            print("Your email doesn't look like an email: '{}'".format(arg))
            print("""Please check it and, if you are sure that your email is right,
  but prober still refuses it  -- contact the developer.""")
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
            print(err_fmt("maximum sequence length must be a positive interger number!"))
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-d", "--indir"):
        if not os.path.isdir(arg):
            print(err_fmt("directory '{}' does not exist!".format(arg)))
            platf_depend_exit(1)
        # end if
        
        indir_path = os.path.abspath(arg)

        # Add all fastq and fasta files from '-d' directory to fq_fa_list
        fq_fa_list.extend(list( filter(is_fq_or_fa, glob("{}{}*".format(indir_path, os.sep))) ))

    elif opt in ("-a", "--algorithm"):
        if not arg in ("0", "1", "2"):
            print(err_fmt("invalid value specified by '-a' option!"))
            print("Available values: 0 for megaBlast, 1 for discoMegablast, 2 for blastn")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end if
        blast_algorithm = ("megaBlast", "discoMegablast", "blastn")[int(arg)]
    # end if
# end for


# If no FASTQ or FASTA file have been found
if len(fq_fa_list) == 0:
    # If input directory was specified -- exit
    if not indir_path is None:
        print(err_fmt("""no input FASTQ or FASTA files specified
  or there is no FASTQ and FASTA files in the input directory.\n"""))
        platf_depend_exit(1)
    
    # If input directory was not specified -- look for FASTQ files in working directory
    else:
        fq_fa_list = list(filter( is_fq_or_fa, glob("{}{}*".format(os.getcwd(), os.sep)) ))

        # If there are nothing to process -- just show help message
        if len(fq_fa_list) == 0:
            print("\nbarapost-prober.py (Version {})\n".format(__version__))
            print("Usage:")
            print("  barapost-prober.py one.fastq.gz another.fasta [...] [OPTIONS]")
            print("For more detailed description, run:")
            print("  barapost-prober.py -h\n")
            platf_depend_exit(0)
        else:
            # Ask if a user wants to proceed or he/she ran it occasionally and wants just help message
            print("\n {} fasta and/or fastq files are found in working directory.\n".format(len(fq_fa_list)))
            error = True
            while error:
                reply = input("""Press ENTER to process them
  or enter 'h' to just see help message:>> """)
                if reply == "":
                    error = False
                    pass
                elif reply == 'h':
                    error = False
                    print('\n' + '-'*15)
                    print("  barapost-prober.py (Version {})\n".format(__version__))
                    print("Usage:")
                    print("  barapost-prober.py one.fastq.gz another.fasta [...] [OPTIONS]")
                    print("For more detailed description, run:")
                    print("  barapost-prober.py -h\n")
                    platf_depend_exit(0)
                else:
                    print("Invalid reply: {}\n".format(reply))
                # end if
            # end while
        # end if
    # end if
# end if

# Check if there are duplicated basenames in input files:
for path in fq_fa_list:
    bname = os.path.basename(path)
    same_bnames = tuple(filter(lambda f: os.path.basename(f) == bname, fq_fa_list))
    if len(same_bnames) != 1:
        print(err_fmt("input files must have different names"))
        print("List of files having same name:")
        for path in same_bnames:
            print("'{}'".format(path))
        # end for
        platf_depend_exit(1)
    # end if
# end for
del bname, same_bnames


# Sort list of files that will be processed -- process them in alphabetical order.
fq_fa_list.sort()
print() # just print new line

from src.filesystem import create_result_directory
from src.filesystem import is_fastq, is_fasta
from src.platform import get_logfile_path
from src.check_connection import check_connection

# If there are any fastq input files -- import fastq record generator
if any(filter(is_fastq, fq_fa_list)):
    from src.fastq import fastq_packets
# end if

# If there are any fasta input files -- import fasta record generator
if any(filter(is_fasta, fq_fa_list)):
    from src.fasta import fasta_packets
# end if

if packet_mode == 0:
    # Make packet size consistent with probing batch size
    packet_size = min(packet_size, probing_batch_size)
# end if

# Create output directory
if not os.path.isdir(outdir_path):
    try:
        os.makedirs(outdir_path)
    except OSError as oserr:
        print(err_fmt("unable to create result directory '{}'".format(outdir_path)))
        print( str(oserr) )
        platf_depend_exit(1)
    # end try
# end if

#                       |===== Proceed =====|

check_connection("https://blast.ncbi.nlm.nih.gov")

logfile_path = get_logfile_path("prober", outdir_path)

printl(logfile_path, "\n |=== barapost-prober.py (version {}) ===|\n".format(__version__))
printl(logfile_path, get_full_time() + "- Start working\n")

from src.prober_modules.prober_spec import look_around
from src.prober_modules.networking import verify_taxids
from src.prober_modules.kernel import submit, retrieve_ready_job

# Make sure that TaxIDs specified by user actually exist
organisms = verify_taxids(taxid_list, logfile_path)

# Print information about the run
printl(logfile_path, " - Output directory: '{}';".format(outdir_path))
printl(logfile_path, " - Logging to '{}'".format(logfile_path))
if user_email != "":
    printl(logfile_path, " - Your email: {}".format(user_email))
# end if
printl(logfile_path, " - Probing batch size: {} sequences;".format("all" if send_all else probing_batch_size))

mode_comment = "number of sequences" if packet_mode == 0 else "sum of sequences' lengths"
printl(logfile_path, " - Packet forming mode: {} ({});".format(packet_mode,mode_comment))
del mode_comment

if packet_mode == 0:
    tmp_str = "sequences"
else:
    tmp_str = "base pairs"
# end if
printl(logfile_path, " - Packet size: {} {};".format(packet_size, tmp_str))
del tmp_str

if max_seq_len < float("inf"):
    printl(logfile_path, " - Maximum length of a sequence to submit: {} bp;".format(max_seq_len))
# end if
printl(logfile_path, " - BLAST algorithm: {};".format(blast_algorithm))
printl(logfile_path, " - Database: nr/nt;")
if len(organisms) > 0:
    for db_slice in organisms:
        printl(logfile_path, "   {};".format(db_slice))
    # end for
# end if

s_letter = '' if len(fq_fa_list) == 1 else 's'
printl(logfile_path, "\n {} file{} will be processed.".format( len(fq_fa_list), s_letter))
with open(logfile_path, 'a') as logfile: # write paths to all input files to log file
    logfile.write("Here they are:\n")
    for i, path in enumerate(fq_fa_list):
        logfile.write("    {}. '{}'\n".format(i+1, path))
    # end for
# end with

printl(logfile_path, '-'*30)


import src.legacy_taxonomy_handling as legacy_taxonomy_handling

# Configure path to file, which will contain info about what GenBank records to download
acc_fpath = os.path.join(outdir_path, "hits_to_download.tsv")

# Configure path to taxonomy directory
taxonomy_dir = os.path.join(outdir_path, "taxonomy")
if not os.path.isdir(taxonomy_dir):
    os.makedirs(taxonomy_dir)
# end if
taxonomy_path = os.path.join(taxonomy_dir, "taxonomy.tsv")

# Check if there is legacy taxonomy file and, if so, reformat it to new (TSV) format
legacy_taxonomy_handling.check_deprecated_taxonomy(outdir_path, logfile_path)


# Dictionary of accessions and record names.
# Accessions are keys, values are tuples: (record_name, number_of_hits).
acc_dict = dict()

# Counter of sequences processed during current run
seqs_processed = [0] # it should be mutable

# Counter of sequences processed concerning possible previous run(s)
glob_seqs_processed = 0

# Varible for stopping execution when probing batch is processed completely.
stop = False

# Further:
# 1. 'previous_data' is a dict of the following structure:
# {
#     "RID": saved_RID <str>,
#     "packet_size_save": saved packet size <int>,
#     "packet_size_mode": saved packet mode <int>,
#     "tsv_respath": path_to_tsv_file_from_previous_run <str>,
#     "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file <int>,
#     "tmp_fpath": path_to_pemporary_file <str>,
#     "decr_pb": valuse decreasing size of probing batch (see below, where this variable is used) <int>
# }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. 'packet' is a dict of the following structure:
#    {
#        "fasta": FASTA_data_containing_query_sequences (str),
#        "qual": dictionary {seq_id: quality}
#    }
#    quality is '-' for fasta files

#                   |===== Kernel loop =====|

# Iterate over input FASTQ and FASTA files:
for i, fq_fa_path in enumerate(fq_fa_list):

    # Create the result directory with the name of FASTQ of FASTA file being processed:
    new_dpath = create_result_directory(fq_fa_path, outdir_path, logfile_path)

    # "hname" means human readable name (i.e. without file path and extention)
    infile_hname = os.path.basename(fq_fa_path)
    infile_hname = re_search(r"(.+)\.(m)?f(ast)?(a|q)(\.gz)?$", infile_hname).group(1) # remove extention

    # Look around and check if there are results of previous runs of this script
    # If 'look_around' returns None -- there is no data from previous run
    previous_data = look_around(outdir_path, new_dpath, fq_fa_path,
        blast_algorithm, acc_dict, probing_batch_size, logfile_path)

    if previous_data is None: # If there is no data from previous run
        num_done_seqs = 0 # no sequences were processed
        tsv_res_path = "{}.tsv".format(os.path.join(new_dpath,
            "classification")) # form path of result tsv file
        tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath,
            infile_hname), blast_algorithm) # form path to temporary file
        saved_RID = None
        saved_packet_size = None
        saved_packet_mode = None
    else: # if there is data from previous run
        num_done_seqs = previous_data["n_done_reads"] # get number of successfully processed sequences
        tsv_res_path = previous_data["tsv_respath"] # result tsv file should be the same as during previous run
        tmp_fpath = previous_data["tmp_fpath"] # temporary file should be the same as during previous run
        saved_RID = previous_data["RID"] # having this RID we can try to get response for last request without resending
        saved_packet_size = previous_data["packet_size_save"]
        saved_packet_mode = previous_data["packet_mode_save"]
        # Let's assume that a user won't modify his/her brobing_batch size between erroneous runs:
        #   subtract num_done_reads if probing_batch_size > num_done_reads.
        probing_batch_size -= previous_data["decr_pb"]
    # end if

    # Take previous run into account
    glob_seqs_processed += num_done_seqs

    if send_all or packet_mode == 1:
        out_of_n = {"msg": "", "npacks": None}
    else:
        # Calculate total number of packets meant to be sent from current FASTA file
        packs_at_all = probing_batch_size // packet_size

        if probing_batch_size % packet_size > 0: # and this is ceiling
            packs_at_all += 1
        # end if

        # Pretty string to print "Request 3 out of 7" if not 'send_all'
        out_of_n = {"msg": " out of {}".format(packs_at_all), "npacks": packs_at_all}
    # enf if

    # Ordinal number of packet meant to be sent (will incrementing after every successful submission)
    pack_to_send = [1] # it should be mutable

    # Choose appropriate record generator:
    packet_generator = fastq_packets if is_fastq(fq_fa_path) else fasta_packets

    # Iterate over packets in current file
    for packet in packet_generator(fq_fa_path, packet_size, num_done_seqs, packet_mode,
        saved_packet_size, saved_packet_mode, max_seq_len, probing_batch_size):

        # Assumption that we need to submit current packet (that we cannot just request for results)
        send = True

        # If current packet has been already send, we can try just to request for results
        if not saved_RID is None:
            send = retrieve_ready_job(saved_RID, packet, packet_size, packet_mode, pack_to_send, seqs_processed,
                fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
                blast_algorithm, user_email, organisms,
                acc_dict, out_of_n)
            saved_RID = None
        # end if

        # Submit current packet to BLAST server
        if send:
            submit(packet, packet_size, packet_mode, pack_to_send, seqs_processed,
                fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
                blast_algorithm, user_email, organisms,
                acc_dict, out_of_n)
        # end if

        if not send_all and seqs_processed[0] >= probing_batch_size: # probing batch is processed -- finish work
            stop = True
            break
        # end if
    # end for
    if stop:
        break
    # end if
# end for
#            |===== End of the kernel loop =====|


# Print some summary:
glob_seqs_processed += seqs_processed[0]
str_about_prev_runs = ", including previous run(s)" if glob_seqs_processed > seqs_processed[0] else ""

printl(logfile_path, '-'*20+'\n')
space_sep_num = "{:,}".format(glob_seqs_processed).replace(',', ' ')
printl(logfile_path, " {} sequences have been processed{}\n".format(space_sep_num,
    str_about_prev_runs))

printl(logfile_path, "Here are Genbank records that can be used for further annotation by 'barapost-local.py'.")
printl(logfile_path, "They are sorted by their occurence in probing batch:")

# Print accessions and record names sorted by occurence
# "-x[1][2]:": minus because we need descending order, [1] -- get tuple of "other information",
#   [2] -- get 3-rd element (occurence)
for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][1]):
    s_letter = "s" if other_info[1] > 1 else ""
    space_sep_num = "{:,}".format(other_info[1]).replace(',', ' ')
    printl(logfile_path, " {} hit{} - {}, '{}'".format(space_sep_num,
        s_letter, acc, other_info[0]))
# end for

# Print number of unkmown sequences, if there are any:
unkn_num = glob_seqs_processed - sum( map(lambda x: x[1], acc_dict.values()) )
if unkn_num > 0:
    s_letter = "s" if unkn_num > 1 else ""
    printl(logfile_path, " {} sequence{} - No significant similarity found".format(unkn_num, s_letter))
# end if

printl(logfile_path, """\nThey are saved in following file:
  '{}'""".format(acc_fpath))
printl(logfile_path, """\nYou can edit this file before running 'barapost-local.py' in order to
  modify list of sequences that will be downloaded from Genbank
  and used as local (i.e. on your local computer) database by 'barapost-local.py'.""")

printl(logfile_path, '\n' + get_full_time() + " - Task is completed\n")
platf_depend_exit(0)