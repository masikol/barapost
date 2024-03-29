# -*- coding: utf-8 -*-
# Module defines finctions that are "miscallaneous" for barapost.

import os
import re
import sys

from xml.etree import ElementTree # for retrieving information from XML BLAST report
import subprocess as sp

from src.printlog import printlog_error, printlog_error_time
from src.platform import platf_depend_exit
from src.filesystem import remove_bad_chars
from src.filesystem import rename_file_verbosely

# Pattern for GenBank accession number
# See https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
#   and maybe https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
GB_ACC_PATTERN = r"([A-Z]{2}_)?([A-Z]{1,2})?[0-9]{5,8}(\.[0-9]+)?"


def look_around(new_dpath, fq_fa_path):
    # Function looks around in order to check if there are results from previous runs of this script.
    #
    # Returns None if there is no result from previous run.
    # If there are results from previous run, returns a dict of the following structure:
    # {
    #     "tsv_respath": path_to_tsv_file_from_previous_run (str),
    #     "n_done_reads": number_of_successfull_requests_from_currenrt_FASTA_file (int),
    # }
    #
    # :param new_dpath: path to current (corresponding to fq_fa_path file) result directory;
    # :type new_dpath: str;
    # :param fq_fa_path: path to current (corresponding to fq_fa_path file) FASTA file;
    # :type fq_fa_path: str;

    # "hname" means human readable name (i.e. without file path and extention)
    fasta_hname = os.path.basename(fq_fa_path) # get rid of absolute path
    fasta_hname = re.search(r"(.*)\.(m)?f(ast)?a", fasta_hname).group(1) # get rid of '.fasta' extention

    # Form path to result file
    tsv_res_fpath = os.path.join(new_dpath, "classification.tsv")

    num_done_reads = 0 # variable to keep number of succeffdully processed sequences

    if os.path.exists(tsv_res_fpath):

        with open(tsv_res_fpath, 'r') as res_file:
            # There can be invalid information in result file
            try:
                lines = res_file.readlines()
                num_done_reads = len(lines) - 1 # the first line is a head
            except OSError as err:
                printlog_error_time("Data in classification file `{}` is broken. Reason:"\
                    .format(tsv_res_fpath))
                printlog_error( str(err) )
                printlog_error("Starting from the beginning.")
                rename_file_verbosely(tsv_res_fpath)
                return None
            # end try
        # end with
    else:
        return None
    # end if

    return {
        "tsv_respath": tsv_res_fpath,
        "n_done_reads": num_done_reads,
    }
# end def look_around


def launch_blastn(packet, blast_algorithm, use_index, queries_tmp_dir, db_path):
    """
    Function launches 'blastn' utility from "BLAST+" toolkit and returns it's response.

    :param pacekt: FASTA data meant to be processend by 'blastn';
    :type packet: str;
    :param blast_algorithm: blastn algorithm to use;
    :type blast_algorithm: str;
    :param use_index: logic value inddicating whether to use index;
    :type use_index: bool:
    :param queries_tmp_dir: path to directory with query files;
    :type queries_tmp_dir: str:
    :param db_path: path to database;
    :type db_path: str:
    """

    # PID of current process won't change, so we can use it to mark query files.
    # 'paket's are too large to pass them to 'subprocess.Popen' as stdin,
    #    therefore we need to use these query files.
    query_path = os.path.join(queries_tmp_dir, "query{}_tmp.fasta".format(os.getpid()))

    with open(query_path, 'w') as query_file:
        query_file.write(packet)
    # end with

    # Configure command line
    blast_cmd = "blastn -query {} -db {} -outfmt 5 -task {} -max_target_seqs 10 -max_hsps 1 -use_index {}"\
        .format(query_path, db_path, blast_algorithm, use_index)

    pipe = sp.Popen(blast_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        printlog_error_time("Error occured while aligning a sequence against local database")
        printlog_error(stdout_stderr[1].decode("utf-8"))
        platf_depend_exit(pipe.returncode)
    # end if

    return stdout_stderr[0].decode("utf-8")
# end def launch_blastn


def parse_align_results_xml(xml_text, qual_dict):
    # Function parses BLAST xml response and returns tsv lines containing gathered information:
    #     1. Query name.
    #     2. Hit name formatted by 'format_taxonomy_name()' function.
    #     3. Hit accession.
    #     4. Length of query sequence.
    #     5. Length of alignment.
    #     6. Percent of identity.
    #     7. Percent of gaps.
    #     8. E-value.
    #     9. Average Phred33 quality of a read (if source file is FASTQ).
    #     10. Read accuracy (%) (if source file is FASTQ).
    #
    # :param xml_text: XML text with results of alignment;
    # :type xml_text: str;
    # :param qual_dict: dict, which maps sequence IDs to their quality;
    # :type qual_dict: dict<str: float>;
    #
    # Returns list<str>.

    result_tsv_lines = list()

    # /=== Parse BLAST XML response ===/

    root = ElementTree.fromstring(xml_text) # get tree instance

    # Iterate over "Iteration" and "Iteration_hits" nodes
    for iter_elem, iter_hit in zip(root.iter("Iteration"), root.iter("Iteration_hits")):

        # "Iteration" node contains query name information
        query_name = iter_elem.find("Iteration_query-def").text
        query_len = iter_elem.find("Iteration_query-len").text

        avg_quality = qual_dict[query_name]
        if avg_quality != '-':
            miscall_prop = round(10**(avg_quality/-10), 3)
            accuracy = round( 100*(1 - miscall_prop), 2 ) # expected percent of correctly called bases
        else:
            # If FASTA file is processing, print dashed in quality columns
            avg_quality = "-"
            accuracy = "-" # expected percent of correctly called bases
        # end if

        # Check if there are any hits
        chck_h = iter_hit.find("Hit")

        if chck_h is None:
            # If there is no hit for current sequence
            result_tsv_lines.append('\t'.join( (query_name, "No significant similarity found", "-", query_len,
                "-", "-", "-", "-", str(avg_quality), str(accuracy)) ))
        else:
            # If there are any hits, node "Iteration_hits" contains at least one "Hit" child
            # Get first-best bitscore and iterato over hits that have the save (i.e. the highest bitscore):
            top_bitscore = next(chck_h.find("Hit_hsps").iter("Hsp")).find("Hsp_bit-score").text

            annotations = list()
            hit_accs = list()

            for hit in iter_hit:

                # Find the first HSP (we need only the first one)
                hsp = next(hit.find("Hit_hsps").iter("Hsp"))

                if hsp.find("Hsp_bit-score").text != top_bitscore:
                    break
                # end if

                curr_acc = sys.intern(hit.find("Hit_accession").text) # get hit accession
                hit_accs.append( curr_acc )

                # Get full hit name (e.g. "Erwinia amylovora strain S59/5, complete genome")
                hit_def = remove_bad_chars(hit.find("Hit_def").text)
                annotations.append(hit_def)

                align_len = hsp.find("Hsp_align-len").text.strip()
                pident = hsp.find("Hsp_identity").text # get number of matched nucleotides
                gaps = hsp.find("Hsp_gaps").text # get number of gaps

                evalue = hsp.find("Hsp_evalue").text # get e-value
            # end for

            # Divide annotations and accessions with '&&'
            annotations = '&&'.join(annotations)
            hit_accs = '&&'.join(hit_accs)

            # Append new tsv line containing recently collected information
            result_tsv_lines.append( '\t'.join( (query_name, annotations, hit_accs, query_len,
                align_len, pident, gaps, evalue, str(avg_quality), str(accuracy)) ))
        # end if
    # end for

    return result_tsv_lines
# end def parse_align_results_xml


def configure_acc_dict(acc_fpath, your_own_fasta_lst, accs_to_download):
    # Fucntion configures accession dictionary according to accession file generated by 'barapost-prober.py':
    #    keys are accessions, values are tuples of the following format:
    #     (<sequence_name_aka_definition>).
    #
    # :param acc_fpath: path to accession file generated by 'barapost-prober.py';
    # :type acc_fpath: str;
    # :param your_own_fasta_lst: list of paths to user's fasta files;
    # :type your_own_fasta_lst: list<str>;
    #
    # Returns accession dictionary described above.

    acc_dict = dict()

    # if database will be created only from 'your own' FASTA files -- return empty dict
    if not acc_fpath is None:

        with open(acc_fpath, 'r') as acc_file:
            lines = acc_file.readlines()

            for line_idx, line in enumerate(lines):
                line = line.strip()
                # Ignore ampty lines, commented lines and head of the table:
                if line != "" and not line.startswith('#') and not line.startswith("ACCESSION"):

                    line_splt = line.split('\t')
                    acc = sys.intern(line_splt[0].partition('.')[0])

                    if not re.match(GB_ACC_PATTERN, acc) is None:
                        # If we encounter GenBank accession number
                        try:
                            if len(line_splt) == 1: # just accession
                                name = "No definition of the sequence provided"
                            else:
                                name = line_splt[1]
                            # end if
                            acc_dict[acc] = name
                        except IndexError as err:
                            printlog_error_time("Error: invalid data in file `{}`!".format(acc_fpath))
                            printlog_error("Here is that invalid line:\n  `{}`".format(line))
                            printlog_error(str(err))
                            platf_depend_exit(1)
                        # end try
                    else:
                        # It it's not a GenBank accession number,
                        #   probably it is a path to reference file.
                        if os.path.exists(line):
                            your_own_fasta_lst.append(line)
                        else:
                            printlog_error_time("Error in file `{}`.".format(acc_fpath))
                            printlog_error("Line #{} looks like path to reference file, but this file does not exist."\
                                .format(line_idx+1))
                            printlog_error("Here is this invalid line:\n  `{}`".format(line))
                            platf_depend_exit(1)
                        # end if
                    # end if
                # end if
            # end for
        # end with
    # end if

    if len(your_own_fasta_lst) == 0 and len(acc_dict) == 0 and len(accs_to_download) == 0:
        printlog_error_time("Error: no accession information found in file `{}`".format(acc_fpath))
        platf_depend_exit(1)
    # end if

    return acc_dict
# end def configure_acc_dict
