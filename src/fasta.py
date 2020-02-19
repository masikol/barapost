# -*- coding: utf-8 -*-
# This module defines function-generator that yields fasta-formatted records from fasta file
#   of certain ('packet_size') size.

from src.prune_seqs import prune_seqs
from src.fmt_readID import fmt_read_id
from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped

def pass_processed_seqs(fasta_file, num_done_seqs, fmt_func):
    """
    Function passes sequences that have been already processed.

    :param fasta_file: FASTA file instalce;
    :type fasta_file: str;
    :param num_done_seqs: amount of sequences that have been already processed;
    :type num_done_seqs: int;
    :param fmt_func: function from 'FORMATTING_FUNCS' tuple;
    """

    if num_done_seqs == 0:
        return None # no sequences have been processed
    else:
        i = 0
        while i <= num_done_seqs:

            line = fmt_func(fasta_file.readline())
            if line == "":
                return ""
            # end if
            if line.startswith('>'):
                line = fmt_read_id(line)
                next_id_line = line
                i += 1
            # end if
        # end while

        # return ID of seqeunce, which will be included in the next packet first
        return next_id_line
    # end if
# end def pass_processed_seqs


def fasta_packets(fasta, packet_size, num_done_seqs, max_seq_len=None):
    """
    Generator yields fasta-formattedpackets of records from fasta file.
    This function passes 'num_done_seqs' sequences (i.e. they will not be processed)
        to 'pass_processed_files'.

    :param fasta: path to fasta file;
    :type fasta: str;
    :param packet_size: number of sequences to align in one request ('blastn' launching);
    :type packet_size: int;
    :param reads_at_all: number of sequences in current file;
    :type reads_at_all: int;
    :param num_done_seqs: number of sequnces in current file that have been already processed;
    :type num_doce_reads: int;
    :param max_seq_len: maximum length of a sequence proessed;
    :type max_seq_len: int (None if pruning is disabled);
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fasta) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fasta) ]

    with how_to_open(fasta) as fasta_file:
        # Next line retrieving is implemented as simple line-from-file reading.
        get_next_line = lambda: fmt_func(fasta_file.readline())

        # Variable that contains ID of next sequence in current FASTA file.
        # If no or all sequences in current FASTA file have been already processed, this variable is None.
        # There is no way to count sequences in multi-FASTA file, accept of counting sequence IDs.
        # Therefore 'next_id_line' should be saved in memory just after moment when packet is formed.
        next_id_line = pass_processed_seqs(fasta_file, num_done_seqs, fmt_func)

        if next_id_line == "":
            yield {"fasta": "", "names": list()}
        # end if

        packet = ""

        # We are resuming, nucleotide sequence will be saved in 'line' variable here:
        line = get_next_line()
        if line.startswith('>'):
            line = fmt_read_id(line) # format sequence ID
        # end if

        # If some sequences have been passed, this if-statement will be executed.
        # New packet should start with sequence ID line.
        if not next_id_line is None:
            packet += next_id_line+'\n'
        # end if
        packet += line+'\n' # add recently read line

        qual_dict = dict() # {<seq_id>: '-'}, as soon as fasta file is being processed
        eof = False

        while not eof: # till the end of file

            i = 0 # variable for counting sequences within packet

            while i < packet_size:

                line = get_next_line()
                if line.startswith('>'):
                    line = fmt_read_id(line)
                    i += 1
                # end if
                
                if line == "": # if end of file (data) is reached
                    break
                # end if
                packet += line + '\n' # add line to packet
            # end while

            if line != "":
                next_id_line = packet.splitlines()[-1] # save sequence ID next packet will start with
                packet = '\n'.join(packet.splitlines()[:-1]) # exclude 'next_id_line' from packet
            else:
                eof = True
                next_id_line = None
            # end if

            # Get list of sequence IDs:
            names = filter(lambda l: True if l.startswith('>') else False, packet.splitlines())
            names = map(lambda l: l.replace('>', ''), names)

            for name in names:
                qual_dict[name] = '-' # there is no quality info in fasta files
            # end for

            if not max_seq_len is None: # prune sequences
                packet = prune_seqs(packet, 'l', max_seq_len)
            # end if

            if packet != "":
                yield {"fasta": packet, "qual": qual_dict}
                # Reset packet
                packet = next_id_line+'\n'
                qual_dict.clear()
            else:
                raise StopIteration
            # end if

        # end while
    # end with
# end def fasta_packets