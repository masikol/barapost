# -*- coding: utf-8 -*-

from src.prune_seqs import prune_seqs
from src.fmt_readID import fmt_read_id

def fasta_packets(data, packet_size, num_done_reads, max_seq_len=None):
    """
    Function (actually, generator) retrieves 'packet_size' records from 'fasta'
        no matter whether is it path to FASTA file of actual FASTA data of 'str' type.
    This function will pass 'num_done_reads' sequences (i.e. they will not be processed)
        by calling 'pass_processed_files'.

    :param data: path to FASTA file OR actual FASTA data of 'str' type;
    :type data: str;
    :param packet_size: number of sequences to align in one 'blastn' launching;
    :type packet_size: int;
    :param reads_at_all: number of sequences in current file;
    :type reads_at_all: int;
    :param num_done_reads: number of sequnces in current file that have been already processed;
    :type num_doce_reads: int;
    """

    fasta_lines = data.splitlines()
    del data # let interpreter get rid of this large string -- we do not need it any more
    line_i = 0 # variable for line counting

    def get_next_line():
        """
        Function retrieves element from 'fasta_lines' by 'line_i' index
            and increases 'line_i' variable by 1.
        """
        nonlocal line_i
        try:
            line = fasta_lines[line_i]
            line_i += 1
        
        except IndexError:
            # if IndexError is raised -- we have reached the end of data.
            # Simulate returning of empty string, just like io.TextIOWrapper.readline() does
            #    if end of file is reached:
            return ""
        
        else:
            return line
        # end try
    # end def get_next_line

    # Variable that contains id of next sequence in current FASTA file.
    # If no or all sequences in current FASTA file have been already processed, this variable is None.
    # There is no way to count sequences in multi-FASTA file, accept of counting sequence IDs.
    # Therefore 'next_id_line' should be saved in memory after moment when packet is formed.
    next_id_line = None

    packet = ""

    line = get_next_line()
    if line.startswith('>'):
        line = fmt_read_id(line) # prune sequence ID
    # end if

    packet += line+'\n' # add recently read line

    qual_dict = dict()
    eof = False

    # Iterate over packets left to process
    while not eof:

        i = 0 # variable for counting sequenes within packet

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
# end def fasta_packets_from_str