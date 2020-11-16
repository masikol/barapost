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


def fasta_packets(fasta, packet_size, num_done_seqs, packet_mode=0,
    saved_packet_size=None, saved_packet_mode=None,
    max_seq_len=float("inf"), probing_batch_size=float("inf")):
    """
    Generator yields fasta-formattedpackets of records from fasta file.
    This function passes 'num_done_seqs' sequences (i.e. they will not be processed)
        to 'pass_processed_files'.

    :param fasta: path to fasta file;
    :type fasta: str;
    :param packet_size: number of sequences to align in one request ('blastn' launching);
    :type packet_size: int;
    :param num_done_seqs: number of sequnces in current file that have been already processed;
    :type num_done_seqs: int;
    :param packet_mode: packet mode (see -c option);
    :type packet_mode: int;
    :param saved_packet_size: size of last sent packet from tmp file. Necessary for resumption.
      It will be None, if no tmp file was in classification directory;
    :type saved_packet_size: int;
    :param saved_packet_mode: mode used whilst formig the last sent packet from tmp file.
      Necessary for resumption. It will be None, if no tmp file was in classification directory;
    :type saved_packet_mode: int;
    :param max_seq_len: maximum length of a sequence proessed;
    :type max_seq_len: int (float("inf") if pruning is disabled);
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
            yield {"fasta": "", "qual": dict()}
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

        # Here goes check for saved packet size and mode:
        if not saved_packet_size is None:
            wrk_pack_size = saved_packet_size
        else:
            wrk_pack_size = packet_size
        # end if

        if not saved_packet_mode is None:
            wrk_pack_mode = saved_packet_mode
        else:
            wrk_pack_mode = packet_mode
        # end if

        eof = False
        while not eof: # till the end of file

            counter = 0 # variable for counting sequences within packet
            seqlen = 0

            while counter < wrk_pack_size:

                line = get_next_line()
                if line.startswith('>'):
                    line = fmt_read_id(line)
                    if packet_mode == 0:
                        counter += 1
                    else:
                        counter += min(seqlen, max_seq_len)
                        seqlen = 0
                    # end if
                # end if
                
                if line == "": # if end of file (data) is reached
                    break
                # end if

                if not line.startswith('>'):
                    seqlen += len(line.strip())
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

            # {<seq_id>: '-'}, as soon as it is a fasta file
            qual_dict = {name: '-' for name in names}

            if max_seq_len < float("inf"): # prune sequences
                packet = prune_seqs(packet, max_seq_len)
            # end if

            if packet != "":
                yield {"fasta": packet, "qual": qual_dict}

                if packet_mode == 0:
                    probing_batch_size -= wrk_pack_size
                    wrk_pack_size = min(packet_size, probing_batch_size)
                else:
                    probing_batch_size -= len(qual_dict)
                # end if

                # Switch back to standart packet size
                # As Vorotos said, repeated assignment is the best check:
                if wrk_pack_mode != packet_mode:
                    wrk_pack_mode = packet_mode
                # end if

                if not next_id_line is None:
                    packet = next_id_line+'\n'
                # end if
            else:
                return
            # end if
        # end while
    # end with
# end def fasta_packets


def fasta_packets_from_str(data, packet_size):
    """
    Generator retrieves 'packet_size' records from 'fasta'
        no matter whether is it path to FASTA file of actual FASTA data of 'str' type.

    :param data: FASTA-formatted string;
    :type data: str;
    :param packet_size: number of sequences to align in one 'blastn' launching;
    :type packet_size: int;
    """

    fasta_lines = data.splitlines()
    fasta_lines.append("") # append fictive value imitating end of file
    del data # let interpreter get rid of this large string -- we do not need it any more

    qual_dict = dict() # {<seq_id>: '-'}, as soon as fasta file is being processed

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

        # Get list of sequence IDs:
        names = filter(lambda l: True if l.startswith('>') else False, packet.splitlines())
        names = map(lambda l: l.replace('>', ''), names)

        for name in names:
            qual_dict[name] = '-' # there is no quality info in fasta files
        # end for

        # No way to get quality from fasta-formatted string.
        # However, we will have it from 'packet_generator()' launching 
        #   (see function 'process' in src/barapost_local_modules/parallel_single_file.py).
        if packet != "":
            yield {"fasta": packet, "qual": qual_dict}
            # Reset packet
            if not next_id_line is None:
                packet = next_id_line+'\n'
                qual_dict = dict()
            # end if
        else:
            return
        # end if
    # end while
# end def fasta_packets_from_str