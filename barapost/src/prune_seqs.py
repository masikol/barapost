# -*- coding: utf-8 -*-
# This module defines function that prunes sequences passed in fasta-formatted string.


def _prune(seq, value):

    len_seq = len(seq)
    if len_seq > value:
        flank = len_seq - value
        amend = flank % 2
        flank //= 2
        seq = seq[flank + amend : (len_seq - flank)]
    # end if
    return seq
# end def _prune


def prune_seqs(packet, value):
    # Function prunes all sequences in fasta-formatted string.
    # Function prunes sequences from both ends.
    # :param packet: FASTA data;
    # :type packet: str;
    # :param value: length of resulting sequence;
    # :type value: int;

    lines = packet.splitlines()
    packet = ""

    if value < 1:
        raise ValueError("Invalid value of parameter `value`: it must be > 1")
    # end if

    # Get indices of ID lines in fasta-formatted string
    id_line_idxs = list(map(lines.index, filter(lambda x: x.startswith('>'), lines)))
    id_line_idxs.append(len(lines)) # fictive last id line

    for id_curr, id_next in zip(id_line_idxs[:-1], id_line_idxs[1:]):
        # If 'i' lies between two ID lines. it's sequence line
        is_seq = lambda i: i > id_curr and i < id_next
        # Get all sequence lines of current record
        seq_lines = map(lambda j: lines[j], filter(is_seq, range(len(lines))))
        seq = _prune("".join(seq_lines), value) # form a single string and prune it
        packet += lines[id_curr] + '\n' + seq + '\n' # apend pruned record to packet
    # end for

    return packet
# end def prune_seqs