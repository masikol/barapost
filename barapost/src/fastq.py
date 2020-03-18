# -*- coding: utf-8 -*-
# This module defines function-generator that yields fasta-formatted records from fastq file
#   of certain ('packet_size') size.
# Moreover, this module can be used to calculate mean quality of a read.

from math import log
from src.prune_seqs import prune_seqs
from src.fmt_readID import fmt_read_id
from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped

FASTQ_LINES_PER_READ = 4

# Function for getting Q value from Phred33 character:
substr_phred33 = lambda q_symb: ord(q_symb) - 33
# List of propabilities corresponding to indices (index is Q, value is the propability):
q2p_map = [10 ** (-q/10) for q in range(128)] # 127 -- max value of a signed byte
# Function for accessing propabilities by Q:
qual2prop = lambda q: q2p_map[q]
# Function for accessing Q by propability:
prop2qual = lambda p: round(-10 * log(p, 10), 2)


def get_read_avg_qual(qual_str):
    """
    Function calculates mean quality of a single read.

    :param qual_str: read's quality line in Phred33;
    :type qual_str: str;
    """

    quals = map(substr_phred33, qual_str) # get Qs
    err_props = map(qual2prop, quals) # convert Qs to propabilities
    avg_err_prop = sum(err_props) / len(qual_str) # calculate average propability
    return prop2qual(avg_err_prop)
# end def get_read_avg_qual


def form_packet(fastq_file, packet_size, fmt_func, max_seq_len):
    """
    Function reads lines from 'fastq_file' and composes a packet of 'packet_size' sequences.

    :param fastq_file: file instance from which to read;
    :type fastq_file: _io.TextIOWrapper or gzip.File;
    :param packet_size: number of sequences to retrive from file;
    :type packet_size: int;
    :param fmt_func: formating functio nfrom FORMMATING_FUNCS tuple;
    :param max_seq_len: maximum length of a sequence proessed;
    :type max_seq_len: int (None if pruning is disabled);
    """

    packet = ""
    qual_dict = dict() # {<seq_id>: <read_quality>}
    eof = False

    for i in range(packet_size):

        read_id = fmt_func(fastq_file.readline())

        if read_id == "": # if eof is reached, leave now
            eof = True
            break
        # end if

        read_id = fmt_read_id(read_id)
        seq = fmt_func(fastq_file.readline())
        fastq_file.readline() # pass comment
        avg_qual = get_read_avg_qual( fmt_func(fastq_file.readline()) )

        packet += read_id + '\n' + seq + '\n'
        qual_dict[read_id[1:]] = avg_qual
    # end for

    if not max_seq_len is None: # prune sequences
        packet = prune_seqs(packet, 'l', max_seq_len)
    # end if

    return {"fasta": packet, "qual": qual_dict}, eof

# end def form_packet


def fastq_packets(fastq, packet_size, num_done_seqs, saved_packet_size=None, max_seq_len=None):
    """
    Generator yields fasta-formattedpackets of records from fastq file.
    This function passes 'num_done_seqs' sequences (i.e. they will not be processed)
      to 'pass_processed_files'.

    :param fasta: path to fastq file;
    :type fasta: str;
    :param packet_size: number of sequences to align in one request ('blastn' launching);
    :type packet_size: int;
    :param reads_at_all: number of sequences in current file;
    :type reads_at_all: int;
    :param num_done_seqs: number of sequnces in current file that have been already processed;
    :type num_doce_reads: int;
    :param saved_packet_size: size of last sent packet from tmp file. Necessary for resumption.
      It will be None, if no tmp file was in classification directory;
    :type saved_packet_size: int;
    :param max_seq_len: maximum length of a sequence proessed;
    :type max_seq_len: int (None if pruning is disabled);
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fastq) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq) ]

    with how_to_open(fastq) as fastq_file:

        # Pass reads, which have been already processed:
        for _ in range(int(num_done_seqs * FASTQ_LINES_PER_READ)):
            fastq_file.readline()
        # end for

        # End of file
        eof = False

        # Here goes check for saved packet size:
        if not saved_packet_size is None:
            tmp_pack_size = saved_packet_size
        else:
            tmp_pack_size = packet_size
        # end if

        # Process all remaining sequences with standart packet size:
        while not eof:
            packet, eof = form_packet(fastq_file, tmp_pack_size, fmt_func, max_seq_len)

            if packet["fasta"] == "":
                return
            # end if

            yield packet

            # Switch back to standart packet size
            # As Vorotos said, repeated assignment is the best check:
            tmp_pack_size = packet_size

            if eof:
                return
            # end if
        # end while
    # end with
# end def fastq_packets