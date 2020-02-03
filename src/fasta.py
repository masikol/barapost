# -*- coding: utf-8 -*-

from src.prune_seqs import prune_seqs
from src.fmt_readID import fmt_read_id

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
        return None
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
        return next_id_line
    # end if
# end def pass_processed_seqs


def fasta_packets(fasta, packet_size, num_done_seqs, max_seq_len=None):
    """
    Generator retrieves 'packet_size' records from FASTA file.
    This function will pass 'num_done_seqs' sequences (i.e. they will not be processed)
        by calling 'pass_processed_files'.

    :param fasta: path to FASTA file;
    :type fasta: str;
    :param packet_size: number of sequences to align in one 'blastn' launching;
    :type packet_size: int;
    :param reads_at_all: number of sequences in current file;
    :type reads_at_all: int;
    :param num_done_seqs: number of sequnces in current file that have been already processed;
    :type num_doce_reads: int;
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fasta) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fasta) ]

    with how_to_open(fasta) as fasta_file:
        # Next line etrieving will be performed as simple line-from-file reading.
        get_next_line = lambda: fmt_func(fasta_file.readline())

        # Variable that contains id of next sequence in current FASTA file.
        # If no or all sequences in current FASTA file have been already processed, this variable is None.
        # There is no way to count sequences in multi-FASTA file, accept of counting sequence IDs.
        # Therefore 'next_id_line' should be saved in memory after moment when packet is formed.
        next_id_line = pass_processed_seqs(fasta_file, num_done_seqs, fmt_func)

        if next_id_line == "":
            yield {"fasta": "", "names": list()}
        # end if

        packet = ""

        line = get_next_line()
        if line.startswith('>'):
            line = fmt_read_id(line) # prune sequence ID
        # end if

        # If some sequences have been passed, this if-statement will be executed.
        # New packet should start with sequence ID line.
        if not next_id_line is None:
            packet += next_id_line+'\n'
        # end if
        packet += line+'\n' # add recently read line

        qual_dict = dict()
        eof = False

        # Iterate over packets left to process
        while not eof:

            i = 0 # variable for counting sequenes within packet

            packet_size = min(packet_size, probing_batch_size - seqs_processed)
            
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
                next_id_line = None
            # end if

            # Get list of sequence IDs:
            names = list( filter(lambda l: True if l.startswith('>') else False, packet.splitlines()) )
            names = list( map(lambda l: l.replace('>', ''), names) )

            for name in names:
                qual_dict[name] = '-'
            # end for
            del names # let it go

            if not max_seq_len is None:
                packet = prune_seqs(packet.strip(), 'l', max_seq_len)
            # end if

            yield {"fasta": packet.strip(), "qual": qual_dict}

            # Reset packet
            if not next_id_line is None:
                packet = next_id_line+'\n'
                qual_dict.clear()
            else:
                return
            # end if
        # end while
    # end with
# end def fasta_packets