# -*- coding: utf-8 -*-
from gzip import open as open_as_gzip
from re import search as re_search
from glob import glob
import os
import sys

is_fastq = lambda f: True if not re_search(r".*\.f(ast)?q(\.gz)?$", f) is None else False

is_gzipped = lambda f: True if f.endswith(".gz") else False
OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line,   # format text line
    lambda line: line.decode("utf-8")  # format gzipped line
)


def err_fmt(text):
    """Function for configuring error messages"""
    return "\n  \a!! - ERROR: " + text + '\n'
# end def print_error


def platf_depend_exit(exit_code):
    """
    A function that asks to press ENTER on Windows
        and exits.

    :type exit_code: int;
    """
    if sys.platform.startswith("win"):
        input("Press ENTER to exit:")
    # end if
    exit(exit_code)
# end def platf_depend_exit


# According to
# https://github.com/nanoporetech/ont_h5_validator/blob/master/h5_validator/schemas/multi_read_fast5.yaml
ont_read_signature = r"([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})"

def fmt_read_id(read_id):
    """Function formats ID of read: it keeps part that matches 'ont_read_signature' if there is one
    and the first 'word' (separator is space) of this string otherwise.

    :param read_id: ID of a read;
    :type read_id: str;
    """
    srch_ont_read = re_search(ont_read_signature, read_id)
    if srch_ont_read is None:
        return read_id.partition(' ')[0].replace('>', '').strip()
    else:
        return srch_ont_read.group(1).strip()
# end def fmt_read_id


class SeqLength:
    """
    Class contains the only one field -- length of a sequence.
    This class is necessary because:
        1) length filtering is disabled by default;
        2) it is useful to create the length-comparing interface that will
            always tell that the sequence is long enough to keep it
            if length filtering is disabled;
    """

    def __init__(self, init_len):
        """
        :param init_len: length value to initialize
        :type init_len: int;
        """
        self.value = init_len
    # end def __init__

    def __lt__(self, rcmp_len):
        """
        "Less than" (<) dunder method.
        It will always return False if length filtering is disabled
            (i.e. if the rigth operand is None).
        :param rcmp_len: the rigth operand of '<' operator;
        :type rcmp_len: int or None;
        """
        if not rcmp_len is None:
            return self.value < rcmp_len
        else:
            return False
        # end if
    # end def __lt__

# end class seq_length


def get_curr_res_dpath(fq_fa_path, tax_annot_res_dir):
    """
    Function configures and returns the path to result directory for particular FASTA or FASTQ file.
    Result directory is nested in 'tax_annot_res_dir' and is named according to name of FASTA/FASTQ file.
    E.g. if file 'some_reads.fastq' is processing, it's result directory will be named 'some_reads'.

    :param fq_fa_path: path to FASTA/FASTQ file that is procesing;
    :type fq_fa_path: str;
    :param tax_annot_res_dir: path to directory with results of 'prober.py';
    :type tax_annot_res_dir: str;

    Returns path to result directory that was recenly created of 'str' type.
    """

    # dpath means "directory path"
    new_dpath = os.path.join(tax_annot_res_dir, os.path.basename(fq_fa_path)) # get rid of absolute path
    new_dpath = re_search(r"(.*)\.(m)?f(ast)?(a|q)", new_dpath).group(1) # get rid of extention

    return new_dpath
# end def get_curr_res_dir


def get_res_tsv_fpath(new_dpath):
    """
    Function returns current TSV file. Sorting will be performad according to this file.

    :param new_dpath: current result directory;
    :type new_dpath: str;
    """

    brpst_resfile_patt = r".*classification\.tsv$"

    is_similar_to_tsv_res = lambda f: True if not re_search(brpst_resfile_patt, f) is None else False

    if not os.path.exists(new_dpath):
        printl(err_fmt("directory '{}' does not exist!".format(new_dpath)))
        printl(""" Please make sure you have performed taxonomic annotation of the following file
    '{}...'
    with 'prober.py' 'barapost.py'""".format(os.path.basename(new_dpath)))
        printl("""Also this error might occur if you forget to specify result directory
    generated by 'prober.py' with '-r' option.""")
        platf_depend_exit(0)
    # end if

    # Recent file will be the first in sorted list
    tsv_res_fpath = list( filter(is_similar_to_tsv_res, sorted(os.listdir(new_dpath))) )[0]

    return os.path.join(new_dpath, tsv_res_fpath)
# end def get_res_tsv_fpath


# There is an accession number in the beginning of local FASTA file
local_name_hit_patt = r"OWN_SEQ_[0-9]+_"
# Pattern that will match ID of seqeunce in FASTA file generated by SPAdes
spades_patt = r"(NODE)_([0-9]+)"
# Pattern that will match ID of seqeunce in FASTA file generated by a5
a5_patt = r"(scaffold)_([0-9]+)"
# Pattern that will match file path in sequence ID
path_patt = r"\(_(.+)_\)"


def format_taxonomy_name(hit_name, sens):
    """
    Function formats taxonomy name according to chosen sensibiliry of sorting.
    :param hit_name: full_fit_name_of_the_subject_sequence;
    :type hit_name: str;
    :param sens: sensibility returned by 'get_classif_sensibility()' function.
        It's value can be one of the following strings: "genus", "sprcies", "strain";
    :type sens: str;
    Returns formatted hit name of 'str' type;
    """

    hit_names = hit_name.split('&&')
    bricks = list()

    for modif_hit_name in hit_names:
        # This string can be edited in this funtion, original hei name will be kept intact
        modif_hit_name = modif_hit_name.strip()

        # If there is no hit -- we are sure what to do!
        if modif_hit_name == "No significant similarity found":
            return "unknown"
        # end if

        # Check if hit is a sequence from SPAdes or a5 assembly:
        spades_match_obj = re_search(spades_patt, modif_hit_name)
        a5_match_obj = re_search(a5_patt, modif_hit_name)

        for match_obj in (spades_match_obj, a5_match_obj):

            # If hit is a sequence from SPAdes or a5 assembly
            if not match_obj is None:

                # Find path to file with assembly:
                try:
                    assm_path = re_search(path_patt, modif_hit_name).group(1)
                except AttributeError:
                    assm_path = None
                # end

                node_or_scaff = match_obj.group(1) # get word "NODE" or "scaffold"
                node_scaff_num = match_obj.group(2) # get it's number

                # SPAdes generate "NODEs"
                if node_or_scaff == "NODE":
                    assmblr_name = "SPAdes"
                # a5 generates "scaffolds"
                elif node_or_scaff == "scaffold":
                    assmblr_name = "a5"
                # There cannot be enything else
                else:
                    print(err_fmt("signature of sequence ID from assembly not recognized: '{}'".format(hit_name)))
                    platf_depend_exit(1)
                # end if

                # Include file path to sorted file name
                # Replace path separetor with underscore in order not to held a bacchanalia in file system.
                if assm_path is not None:
                    if sens == "genus":
                        # Return only path and "NODE" in case of SPAdes and "scaffold" in case of a5
                        bricks.append('_'.join( (assmblr_name, "assembly", assm_path.replace(os.sep, '_'), node_or_scaff) ))
                    else:
                        # Return path and "NODE_<N>" in case of SPAdes and "scaffold_<N>" in case of a5
                        bricks.append('_'.join( (assmblr_name, "assembly", assm_path.replace(os.sep, '_'), node_or_scaff,
                            node_scaff_num) ))
                    # end if
                else:
                    if sens == "genus":
                        # Return only "NODE" in case of SPAdes and "scaffold" in case of a5
                        bricks.append(assmblr_name + "_assembly_" + node_or_scaff)
                    else:
                        # Return "NODE_<N>" in case of SPAdes and "scaffold_<N>" in case of a5
                        bricks.append('_'.join( (assmblr_name + "assembly", node_or_scaff, node_scaff_num )))
                    # end if
                # end if
            # end if
        # end for

        if not ' ' in modif_hit_name:
            # We have lineage
            taxa_splitnames = modif_hit_name.split(';')
            if not re_search(r";[a-z]", modif_hit_name) is None:
                if sens == "genus":
                    bricks.append(taxa_splitnames[-2])
                else:
                    bricks.append(taxa_splitnames[-2] + '_' + taxa_splitnames[-1])
                # end if
            else:
                try:
                    if sens == "genus":
                        bricks.append(taxa_splitnames[5])
                    else:
                        raise IndexError
                    # end if
                except IndexError:
                    bricks.append("unknown")
                # end try
            # end if
        else:

            #                                      Genus    species                   strain name and anything after it
            hit_name_patt = r"(PREDICTED)?(:)?(_)?[A-Z][\.a-z]+_[a-z]*(sp\.)?(phage)?_(strain_)?.+$"

            # If structure of hit name is strange
            if re_search(hit_name_patt, modif_hit_name) is None:
                bricks.append(modif_hit_name.strip().replace(' ', '_'))   # return full name
            # end if

            taxa_name = modif_hit_name.partition(',')[0]
            taxa_splitnames = taxa_name.strip().split('_')

            # Sometimes query sequence hits records lke this:
            # XM_009008688, 'PREDICTED: Callithrix jacchus cyclin dependent kinase inhibitor 1C (CDKN1C), mRNA'
            if "PREDICTED" in taxa_splitnames[0].upper():
                taxa_splitnames = taxa_splitnames[1:]
            # end if

            # If hit is a phage sequence
            if taxa_splitnames[1] == "phage":
                # Assumming that the man who sortes by genus or species isn't interested in phage strain name
                bricks.append('_'.join( [taxa_splitnames[0], taxa_splitnames[1]] )) # return "<Host_name> phage"
            # end if

            # E.g. 'Bacterium clone zdt-9n2'
            if "clone" in taxa_splitnames:
                bricks.append(taxa_name.replace(' ', '_'))   # return full name
            # end if

            # If someone has shortened genus name
            if '.' in taxa_splitnames[0] or '.' in taxa_splitnames[1]:
                # 'E. coli'
                if not re_search(r"^[A-Z]\.$", taxa_splitnames[0]) is None:
                    bricks.append('_'.join( [taxa_splitnames[0], taxa_splitnames[1]] )) # return genus and species
                # 'E.coli' (without space)
                elif not re_search(r"^[A-Z]\.[a-z]+$", taxa_splitnames[0]) is None:
                    bricks.append(taxa_splitnames[0]) # 'E.coli' will be returned
            # end if 

            if sens == "genus":
                bricks.append(taxa_splitnames[0]) # return genus
            
            elif sens == "species":
                # if species is not specified
                if taxa_splitnames[1] == "sp.":
                    return taxa_name.replace(' ', '_')   # return full name
                else:
                    bricks.append('_'.join( [taxa_splitnames[0], taxa_splitnames[1]] )) # return genus and species
                # end if
            # end if
        # end if
    # end for

    return "&&".join(bricks)
# end def format_taxonomy_name


def fastq_records(fq_path):
    """
    :param read_file: file instance of FASTQ file to retrieve sequences from;
    :type fasta_file: _io.TextIOWrapper or gzip.GzipFile;
    :param fmt_func: function from 'FORMATTING_FUNCS' tuple;

    Returns dictionary of the following structure:
    {
        "seq_id": ID_of_sequence,
        "seq": sequence_itself,
        "opt_id": the_third_line,
        "qual_line": quality_line
    }
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fq_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_path) ]
    # count number of reads in file (4 lines per read)
    num_reads = int(sum(1 for line in how_to_open(fq_path) if line != '') / 4)

    with how_to_open(fq_path) as fq_file:

        for _ in range(num_reads):
            # Read 4 lines of fastq-record
            yield {
                "seq_id": fmt_func(fq_file.readline()),
                "seq": fmt_func(fq_file.readline()),
                "opt_id": fmt_func(fq_file.readline()),
                "qual_line": fmt_func(fq_file.readline())
            }
        # end for
    # end with
    return
# end def read_fastq_record

def fasta_records(fa_path):
    """
    :param fasta_file: file instance of FASTA file to retrieve sequences from;
    :type fasta_file: _io.TextIOWrapper or gzip.GzipFile;
    :param fmt_func: function from 'FORMATTING_FUNCS' tuple;

    Returns dictionary of the following structure:
    {
        "seq_id": ID_of_sequence,
        "seq": sequence_itself
    }
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fa_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fa_path) ]

    with how_to_open(fa_path) as fa_file:

        line = fmt_func(fa_file.readline())
        seq_id = line
        seq = ""
        line = fmt_func(fa_file.readline())

        while line != "":

            seq += line
            line = fmt_func(fa_file.readline())

            if line.startswith('>') or line == "":

                yield {"seq_id": seq_id, "seq": seq}
                seq_id = line
                seq = ""
                line = fmt_func(fa_file.readline())
            # end if
        # end while
    # end with
    return
# end def read_fasta_record


def configure_resfile_lines(tsv_res_fpath, sens):
    """
    Function returns dictionary, where keys are sequence (i.e. sequences meant to be sorted) IDs,
        and values are corresponding hit names.

    :param tsv_res_fpath: path to current TSV file. Sorting will be performed accorfing to this TSV file;
    :type tsv_res_fpath: str;
    """

    resfile_lines = dict()

    with open(tsv_res_fpath, 'r') as brpst_resfile:

        brpst_resfile.readline() # pass the head of the table
        line = brpst_resfile.readline().strip() # get the first informative line

        while line != "":
            splt = line.split('\t')
            read_name = sys.intern(splt[0])
            hit_name = splt[1]
            try:
                query_len = int(splt[3])  # we will filter by length
            except ValueError as verr:
                printl(err_fmt("query length parsing error"))
                printl( str(verr) )
                printl("Please, contact the developer.")
                platf_depend_exit(1)
            # end try
            try:
                phred33 = float(splt[8]) # we will filter by quality
            except ValueError as verr:
                if splt[8] == '-':
                    # Keep minus as quality if there is no quality information.
                    # Error will not be raised.
                    phred33 = splt[8]
                else:
                    printl(err_fmt("query quality parsing error"))
                    printl( str(verr) )
                    printl("Please, contact the developer.")
                    platf_depend_exit(1)
                # end if
            # end try
            resfile_lines[read_name] = [hit_name, phred33, query_len]
            line = brpst_resfile.readline().strip() # get next line
        # end while
    # end with

    # |===== Format taxonomy names =====|
    for read_name in resfile_lines.keys():
        resfile_lines[read_name][0] = format_taxonomy_name(resfile_lines[read_name][0], sens)
    # end for

    return resfile_lines
# end def configure_resfile_lines


def get_checkstr(fast5_fpath):
    """
    Function returns string that will help fasQA5-sorter to find
        TSV file generated by prober and barapost while provessing FASTQ file
        that in turn is basecalled 'fast5_fpath'-file.
    
    Function first searches for ID given to file by programs like of MinKNOW.
    That is:
        1) sequence of 40 (I've seen 40, maybe there can be other number)
        latin letters in lower case interlaced with numbers;
        2) underscore;
        3) number of file within sequenator run;
    For example: file named "FAK94973_e6f2851ddd414655574208c18f2f51e590bf4b27_0.fast5"
        has checkstring "e6f2851ddd414655574208c18f2f51e590bf4b27_0".
    "FAK94973" is not regarding because it can be pruned by basecaller. For example, Guppy acts in this way.

    If no such distinctive string is found in FAST5 file name
        (file can be renamed by the user after sequensing)
        whole file name (except of the '.fast5' extention) is returned as checksting.

    :param fast5_fpath: path to FAST5 file meant to be processed;
    :type fast5_fpath: str;

    Returns checkstring described above.
    """

    try:
        # I'll lower the 40-character barrier down to 30 just in case.
        filename_payload = re_search(r"([a-zA-Z0-9]{30,}_[0-9]+)", fast5_fpath).group(1)
    except AttributeError:
        return os.path.basename(fast5_fpath).replace(".fast5", "")
    else:
        return filename_payload
    # end try
# end def get_checkstr
