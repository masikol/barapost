# -*- coding: utf-8 -*-
# This packet defines generator that splits fasta-formatted string
#   into smaller "packets" and yields them.

from src.prune_seqs import prune_seqs
from src.fmt_readID import fmt_read_id


def fasta_packets_from_str(data, packet_size, max_seq_len=None):
    """
    Generator retrieves 'packet_size' records from 'fasta'
        no matter whether is it path to FASTA file of actual FASTA data of 'str' type.

    :param data: FASTA-formatted string;
    :type data: str;
    :param packet_size: number of sequences to align in one 'blastn' launching;
    :type packet_size: int;
    :param reads_at_all: number of sequences in current file;
    :type reads_at_all: int;
    :param max_seq_len: maximum length of a sequence proessed;
    :type max_seq_len: int (None if pruning is disabled);
    """

    fasta_lines = data.splitlines()
    fasta_lines.append("") # append fictive value imitating end of file
    del data # let interpreter get rid of this large string -- we do not need it any more

    # Variable for counting lines (it is list in order to circumvent immutability of int type)
    line_i = 0

    # Variable that contains id of next sequence in current FASTA file.
    # If no or all sequences in current FASTA file have been already processed, this variable is None.
    # There is no way to count sequences in multi-FASTA file, accept of counting sequence IDs.
    # Therefore 'next_id_line' should be saved in memory after moment when packet is formed.
    next_id_line = None

    # The first line is always a sequence ID if data is not a file, but a fasta-formatted string.
    #   Because in this case all "done" sequences are already passed by function 'fasta_packets'
    line = fmt_read_id(fasta_lines[line_i])
    line_i += 1

    packet = line+'\n' # add recently read line

    eof = False

    # Iterate over packets left to process
    while not eof:

        i = 0 # variable for counting sequenes within packet

        while i < packet_size:

            line = fasta_lines[line_i]

            if line == "": # if end of data is reached
                eof = True
                break
            # end if

            line_i += 1

            if line.startswith('>'):
                line = fmt_read_id(line)
                i += 1
            # end if
            packet += line + '\n' # add line to packet
        # end while

        if line != "":
            next_id_line = packet.splitlines()[-1] # save sequence ID next packet will start with
            packet = '\n'.join(packet.splitlines()[:-1]) # exclude 'next_id_line' from packet
        else:
            next_id_line = None
        # end if

        if not max_seq_len is None:
            packet = prune_seqs(packet.strip(), 'l', max_seq_len)
        # end if

        # No way to get quality from fasta-formatted string.
        # However, we will have it from 'packet_generator()' launching 
        #   (see function 'process' in src/barapost_modules/parallel_single_file.py).
        if packet != "":
            yield {"fasta": packet}
            # Reset packet
            if not next_id_line is None:
                packet = next_id_line+'\n'
            # end if
        else:
            return
        # end if
    # end while
# end def fasta_packets_from_str