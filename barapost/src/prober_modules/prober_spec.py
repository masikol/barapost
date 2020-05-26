# -*- coding: utf-8 -*-
# Module defines finctions that are "miscallaneous" for prober.

import os
import sys
import shelve
from xml.etree import ElementTree # for retrieving information from XML BLAST report
from re import search as re_search

from src.printlog import printl, err_fmt, println, getwt
from src.platform import platf_depend_exit
from src.filesystem import OPEN_FUNCS, rename_file_verbosely, is_gzipped

from src.lineage import find_lineage


def ask_for_resumption(logfile_path):
    """
    Function asks a user if he/she wants to resume the previous run.

    :param logfile_path: path to logfile;
    :type logfile_path: str;

    Returns True if the decision is to resume, else False
    """

    resume = None

    while resume is None:
        resume = input("""
Would you like to resume the previous run?
   1 -- Resume!
   2 -- Start from the beginning.

Enter a number (1 or 2):>> """)
        # Check if entered value is integer number. If no, give another attempt.
        try:
            resume = int(resume)
            # Check if input number is 1 or 2
            if resume != 1 and resume != 2:
                print("\n   Not a VALID number entered!\a\n" + '~'*20)
                resume = None
            else:
                action = "resume the previous run" if resume == 1 else "start from the beginning"
                printl(logfile_path, "You have chosen to {}.\n".format(action))
            # end if
        except ValueError:
            print("\nNot an integer number entered!\a\n" + '~'*20)
            resume = None
        # end try

    return True if resume == 1 else False
# end def ask_for_resumption


def look_around(outdir_path, new_dpath, infile_path, blast_algorithm, acc_dict, probing_batch_size, logfile_path):
    """
    Function looks around in order to check if there are results from previous run(s) of this script
      in order to resume the previous run.

    Returns None if there is no result from previous run.
    If there are results from previous run, returns a dict of the following structure:
    {
        "RID": saved_RID <str>,
        "packet_size_save": saved packet size <int>,
        "packet_size_mode": saved packet mode <int>,
        "tsv_respath": path_to_tsv_file_from_previous_run <str>,
        "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file <int>,
        "tmp_fpath": path_to_pemporary_file <str>,
        "decr_pb": valuse decreasing size of probing batch (see below, where this variable is defined) <int>
    }
    
    :param outdir_path: path to output directory;
    :type outdir_path: str;
    :param new_dpath: path to current (corresponding to fq_fa_path file) result directory;
    :type new_dpath: str;
    :param infile_path: path to current (corresponding to fq_fa_path file) FASTA file;
    :type infile_path: str;
    :param blast_algorithm: BLASTn algorithm to use.
        This parameter is necessary because it is included in name of result files;
    :param acc_dict: dictionary of accession info of hits;
    :type acc_dict: dict<str: tuple<str, str, int>>;
    :param probing_batch_size: amount of sequences meant to be processed in a single run;
    :type probing_batch_size: str;
    :type blast_algorithm: str;
    :param logfile_path: path to logfile;
    :type logfile_path: str;
    """

    # "hname" means human readable name (i.e. without file path and extention)
    fasta_hname = os.path.basename(infile_path) # get rid of absolute path
    fasta_hname = re_search(r"(.*)\.(m)?f(ast)?a", fasta_hname).group(1) # get rid of '.fasta' extention
    how_to_open = OPEN_FUNCS[ is_gzipped(infile_path) ]

    # Form path to temporary file
    tmp_fpath = "{}_{}_temp.txt".format(os.path.join(new_dpath, fasta_hname), blast_algorithm)
    # Form path to result file
    tsv_res_fpath = os.path.join(new_dpath, "classification.tsv")
    # Form path to file with hits to download
    acc_fpath = os.path.join(outdir_path, "hits_to_download.tsv")

    num_done_seqs = 0 # variable to keep number of successfully processed sequences

    resume = None
    # Check if there are results from previous run.
    if os.path.exists(tsv_res_fpath) or os.path.exists(tmp_fpath):
        printl(logfile_path, '\n' + getwt() + " - A result file from previous run is found in the directory:")
        printl(logfile_path, "   '{}'".format(new_dpath))
        # Allow politely to continue from last successfully sent packet.
        resume = ask_for_resumption(logfile_path)
    # end if

    if resume == False:
        rename_file_verbosely(tsv_res_fpath, logfile_path)
        rename_file_verbosely(tmp_fpath, logfile_path)
        rename_file_verbosely(acc_fpath, logfile_path)
    elif resume == True:
        printl(logfile_path, "Let's try to resume...")

        # Collect information from result file
        if os.path.exists(tsv_res_fpath):
            # There can be invalid information in this file
            try:
                with open(tsv_res_fpath, 'r') as res_file:
                    lines = res_file.readlines()
                    num_done_seqs = len(lines) - 1 # the first line is a head
                    last_line = lines[-1]
                    last_seq_id = last_line.split('\t')[0]
                # end with
                # There must be 10 columns in each row:
                if any(map(lambda l: l.count('\t') != 9, lines)):
                    raise ValueError("There must be 10 colums separated by tabs in file 'classification.tsv'")
                # end if

            except Exception as err:
                printl(logfile_path, "\nData in classification file '{}' not found or broken. Reason:".format(tsv_res_fpath))
                printl(logfile_path, ' ' + str(err))

                # If the reason is known -- print erroneous lines
                if isinstance(err, ValueError):
                    printl(logfile_path, "Here are numbers of improper lines:")
                    for i, line in enumerate(lines):
                        if line.count('\t') != 9:
                            printl(logfile_path, str(i+1) + r": '{}'".format(line))
                        # end if
                    # end for
                    printl(logfile_path)
                # end if

                # Ask a user if he/she wants to start from the beginning or to quit
                error = True
                while error:
                    reply = input("""Press ENTER to start from the beginning
  or enter 'q' to quit:>> """)
                    if reply == "":
                        error = False
                        printl(logfile_path, "You have chosen to start from the beginning.\n")
                        rename_file_verbosely(tsv_res_fpath, logfile_path)
                        rename_file_verbosely(tmp_fpath, logfile_path)
                        rename_file_verbosely(acc_fpath, logfile_path)
                        return None
                    elif reply == 'q':
                        platf_depend_exit(0)
                    else:
                        print("! - Invalid reply: '{}'\n".format(reply))
                    # end if
                # end while
            else:
                printl(logfile_path, "Last classified sequence: " + last_seq_id)
                printl(logfile_path, "{} sequences have been already processed".format(num_done_seqs))
            # end try
        # end if

        # Collect information from accession file
        if os.path.exists(acc_fpath):

            # There can be invalid information in this file
            try:
                with open(acc_fpath, 'r') as acc_file:
                    lines = acc_file.readlines()[9:] # omit description and head of the table
                    local_files_filtered = list( filter(lambda x: False if os.path.exists(x) else True, lines) ) # omit file paths
                    for line in local_files_filtered:
                        vals = line.split('\t')
                        acc = sys.intern(vals[0].strip())
                        if len(vals) == 3:
                            acc_dict[acc] = [ vals[1].strip(), int(vals[2].strip()) ]
                        elif len(vals) == 4: # backward compatibility
                            acc_dict[acc] = [ vals[2].strip(), int(vals[3].strip()) ]
                        else:
                            raise IndexError
                        # end if
                    # end for
                # end with

            except Exception as err:
                printl(logfile_path, "\nData in accession file '{}' not found or broken. Reason:".format(acc_fpath))
                printl( logfile_path, ' ' + str(err) )

                # If the reason is known -- print erroneous lines
                if isinstance(err, IndexError):
                    printl(logfile_path, "Here are numbers of improper lines:")
                    for i, line in enumerate(lines):
                        if line.count('\t') != 3:
                            printl(logfile_path, str(i+1) + ": '{}'".format(line))
                        # end if
                    # end for
                    printl(logfile_path)
                # end if

                # Ask a user if he/she wants to start from the beginning or to quit
                error = True
                while error:
                    reply = input("""Press ENTER to start from the beginning
  or enter 'q' to quit:>> """)
                    if reply == "":
                        error = False
                        printl(logfile_path, "You have chosen to start from the beginning.\n")
                        rename_file_verbosely(tsv_res_fpath, logfile_path)
                        rename_file_verbosely(tmp_fpath, logfile_path)
                        rename_file_verbosely(acc_fpath, logfile_path)
                        return None
                    elif reply == 'q':
                        platf_depend_exit(0)
                    else:
                        print("! - Invalid reply: '{}'\n".format(reply))
                    # end if
                # end while
            else:
                printl(logfile_path, "\nHere are Genbank records encountered during previous run:")
                for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][1]):
                    s_letter = "s" if other_info[1] > 1 else ""
                    printl(logfile_path, " {} hit{} - {}, '{}'".format(other_info[1], s_letter, acc, other_info[0]))
                # end for
                print('-'*20)
            # end try
        # end if

        # Get packet size, number of the last sent packet and RID from temp file.
        # There can be invalid information in tmp file of tmp file may not exist
        try:

            with open(tmp_fpath, 'r') as tmp_file:
                temp_lines = tmp_file.readlines()
            # end with

            RID_save = re_search(r"Request_ID: (.+)", temp_lines[0]).group(1).strip()
            packet_size_save = int(re_search(r"Packet_size: ([0-9]*)", temp_lines[1]).group(1).strip())
            packet_mode_save = int(re_search(r"Packet_mode: ([0-9]{1})", temp_lines[2]).group(1).strip())

        except (AttributeError, OSError):

            # There is no need to disturb a user, merely proceed.
            return {
                "RID": None,
                "packet_size_save": None,
                "packet_mode_save": None,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_seqs,
                "tmp_fpath": tmp_fpath,
                "decr_pb": 0
            }
        else:
            # Let's assume that a user won't modify his/her brobing_batch size between erroneous runs:
            #   subtract num_done_reads if probing_batch_size > num_done_reads.
            decr_pb = num_done_seqs if num_done_seqs < probing_batch_size else 0
            # Return data from previous run
            return {
                "RID": RID_save,
                "packet_size_save": packet_size_save,
                "packet_mode_save": packet_mode_save,
                "tsv_respath": tsv_res_fpath,
                "n_done_reads": num_done_seqs,
                "tmp_fpath": tmp_fpath,
                "decr_pb": decr_pb
            }
        # end try
    # end if

    return None
# end def look_around


def parse_align_results_xml(xml_text, qual_dict, acc_dict, logfile_path, taxonomy_path):
    """
    Function parses BLAST xml response and returns tsv lines containing gathered information:
      1. Query name.
      2. Hit name formatted by 'format_taxonomy_name()' function.
      3. Hit accession.
      4. Length of query sequence.
      5. Length of alignment.
      6. Percent of identity.
      7. Percent of gaps.
      8. E-value.
      9. Average quality of a read (if source file is FASTQ).
      10. Read accuracy (%) (if source file is FASTQ).

    :param xml_text: XML text with results of alignment;
    :type xml_text: str;
    :param qual_dict: dict, which maps sequence IDs to their quality;
    :type qual_dict: dict<str: float>;
    :param acc_dict: dictionary comntaining accession data of hits;
    :type acc_dict: dict<str: tuple<str, str, int>>;
    :param logfile_path: path to logfile;
    :type logfile_path: str;
    :param taxonomy_path: path to DBM file with taxonomy;
    :type taxonomy_path: str;

    Returns list<str>.
    """

    result_tsv_lines = list()

    # /=== Parse BLAST XML response ===/

    root = ElementTree.fromstring(xml_text) # get tree instance

    with shelve.open(taxonomy_path, 'c') as tax_file:

        # Iterate through "Iteration" and "Iteration_hits" nodes
        for iter_elem, iter_hit in zip(root.iter("Iteration"), root.iter("Iteration_hits")):
            # "Iteration" node contains query name information
            query_name = sys.intern(iter_elem.find("Iteration_query-def").text)
            query_len = iter_elem.find("Iteration_query-len").text

            avg_quality = qual_dict[query_name]
            if avg_quality != '-':
                miscall_prop = round(10**(avg_quality/-10), 3)
                accuracy = round( 100*(1 - miscall_prop), 2 ) # expected percent of correctly called bases
                qual_info_to_print = "    Average quality of this read is {}, i.e. accuracy is {}%;\n".format(avg_quality,
                    accuracy)
            else:
                # If FASTA file is processing, print dashed in quality columns
                avg_quality = "-"
                accuracy = "-" # expected percent of correctly called bases
                qual_info_to_print = ""
            # end if

            # Check if there are any hits
            chck_h = iter_hit.find("Hit")

            if chck_h is None:
                # If there is no hit for current sequence
                printl(logfile_path, "\n '{}' -- No significant similarity found;\n    Query length - {};".format(query_name, query_len))
                result_tsv_lines.append('\t'.join( (query_name, "No significant similarity found", "-", query_len,
                    "-", "-", "-", "-", str(avg_quality), str(accuracy)) ))
            else:
                # If there are any hits, node "Iteration_hits" contains at least one "Hit" child
                # Get first-best bitscore and iterato over hits that have the save (i.e. the highest bitscore):
                top_bitscore = next(chck_h.find("Hit_hsps").iter("Hsp")).find("Hsp_bit-score").text

                annotations = list()
                hit_accs = list()

                for hit in iter_hit:

                    # Find the first HSP
                    hsp = next(hit.find("Hit_hsps").iter("Hsp"))

                    if hsp.find("Hsp_bit-score").text != top_bitscore:
                        break
                    # end if

                    # Get full hit name (e.g. "Erwinia amylovora strain S59/5, complete genome")
                    hit_def = hit.find("Hit_def").text.replace(' ', '_')
                    annotations.append(hit_def)

                    curr_acc = sys.intern(hit.find("Hit_accession").text)
                    hit_accs.append( curr_acc ) # get hit accession

                    # Get taxonomy
                    find_lineage(curr_acc, hit_def, tax_file,
                        logfile_path)

                    # Update accession dictionary
                    try:
                        acc_dict[curr_acc][1] += 1
                    except KeyError:
                        acc_dict[curr_acc] = [hit_def, 1]
                    # end try

                    align_len = hsp.find("Hsp_align-len").text.strip()
                    pident = hsp.find("Hsp_identity").text # get number of matched nucleotides
                    gaps = hsp.find("Hsp_gaps").text # get number of gaps

                    evalue = hsp.find("Hsp_evalue").text # get e-value
                    pident_ratio = round( float(pident) / int(align_len) * 100, 2)
                    gaps_ratio = round( float(gaps) / int(align_len) * 100, 2)
                # end for

                # Divide annotations and accessions with '&&'
                annotations = '&&'.join(annotations)
                hit_accs = '&&'.join(hit_accs)

                printl(logfile_path, """\n{} - {}
    Query length - {} nt;
    Identity - {}/{} ({}%); Gaps - {}/{} ({}%);""".format(query_name, annotations,
                    query_len, pident, align_len, pident_ratio, gaps, align_len, gaps_ratio))

                # Append new tsv line containing recently collected information
                result_tsv_lines.append( '\t'.join( (query_name, annotations, hit_accs, query_len,
                    align_len, pident, gaps, evalue, str(avg_quality), str(accuracy)) ))

            # end if
            println(logfile_path, qual_info_to_print)
        # end for
    # end with

    return result_tsv_lines
# end def parse_align_results_xml


def write_hits_to_download(acc_dict, acc_file_path):
    """
    Function writes accession data of hits to file 'hits_to_download.tsv'.

    :param acc_dict: dictionary comntaining accession data of hits;
    :type acc_dict: dict<str: tuple<str, str, int>>;
    :param acc_file_path: path to file 'hits_to_download.tsv';
    :type acc_file_path: str;
    """

    # === Write accession information ===

    # Rewrite head:
    with open(acc_file_path, 'w') as acc_file:
        acc_file.write("# Here are accessions and names of Genbank records that can be used for sorting by 'barapost.py'\n")
        acc_file.write("# Values in this file are delimited by tabs.\n")
        acc_file.write("# You are welcome to edit this file by adding,\n")
        acc_file.write("#   removing or muting lines (with adding '#' symbol in it's beginning, just like this description).\n")
        acc_file.write("# Lines muted with '#' won't be noticed by 'barapost.py'.\n")
        acc_file.write("# You can specify your own FASTA files that you want to use as database for 'barapost.py'.\n")
        acc_file.write("# To do it, just write your FASTA file's path to this TSV file in new line.\n\n")
        acc_file.write('\t'.join( ["ACCESSION", "RECORD_NAME", "OCCURRENCE_NUMBER"] ) + '\n')

        # Append updated information:
        for acc, other_info in sorted(acc_dict.items(), key=lambda x: -x[1][1]):
            acc_file.write('\t'.join( (acc, other_info[0], str(other_info[1]))) + '\n')
        # end for
    # end with
# end def write_hits_to_download