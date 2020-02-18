# -*- coding: utf-8 -*-
# This module defines function that prunes sequences passed in fasta-formatted string.

def prune_seqs(packet, mode, value):
    """
    Function prunes all sequences in fasta-formatted string, leaving 5'-half of the sequence.

    :param packet: FASTA data;
    :type packet: str;
    :param mode: available values: 'l' (for length) and 'f' (for fraction).
        Meaning: if mode is 'l', all sequences are pruned from the 5'-end to position 'value' (value > 1),
                 if mode is 'f', all sequences are pruned from the 5'-end to position L*value (0<value<1),
                 where L denotes to length of a sequence.
    :type mode: str;
    :param value: see description if 'mode' parameter above;
    :type value: int or float;
    """

    lines = packet.splitlines()
    packet = ""

    # Define different pruning functions for different modes:
    if mode == 'l':
        if value < 1:
            raise ValueError("Invalid value of parameter 'value': it must be > 1")
        # end if

        def prune(seq, value):
            try:
                return seq[: value]
            except IndexError:
                return seq
            # end try
        # end def prune
    elif mode == 'f':
        if value < 0.0 + 1e-9 or value > 1.0 - 1e-9:
            raise ValueError("Invalid value of parameter 'value': it must be in interval (0, 1)")
        # end if

        def prune(seq, value):
            return seq[: int(len(seq) * value)]
        # end def prune
    else:
        raise ValueError("Invalid 'mode' argument passed to function prune_seqs: '{}'".format(mode))
    # end if

    # Get indices of ID lines in fasta-formatted string
    id_line_idxs = list(map(lines.index, filter(lambda x: x.startswith('>'), lines)))
    id_line_idxs.append(len(lines)) # fictive last id line

    for id_curr, id_next in zip(id_line_idxs[:-1], id_line_idxs[1:]):
        # If 'i' lies between two ID lines. it's sequence line
        is_seq = lambda i: i > id_curr and i < id_next
        # Get all sequence lines of current record
        seq_lines = map(lambda j: lines[j], filter(is_seq, range(len(lines))))
        seq = "".join(seq_lines) # form a single string
        seq = prune(seq, value) # prune it
        packet += lines[id_curr] + '\n' + seq + '\n' # apend pruned record to packet
    # end for

    return packet
# end def prune_seqs