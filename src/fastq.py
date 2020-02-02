# -*- coding: utf-8 -*-

from math import log
from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped
from src.prune_seqs import prune_seqs
from src.fmt_readID import fmt_read_id

FASTQ_LINES_PER_READ = 4

# Function for getting Q value from Phred33 character:
substr_phred33 = lambda q_symb: ord(q_symb) - 33
# List of propabilities corresponding to indices (index is Q, is the propability):
q2p_map = [10 ** (-q/10) for q in range(128)] # 127 -- max value of a signed byte
# Function for accessing propabilities by Q:
qual2prop = lambda q: q2p_map[q]
# Function for accessing Q by propability:
prop2qual = lambda p: round(-10 * log(p, 10), 2)

def get_read_avg_qual(qual_str):

    quals = map(substr_phred33, qual_str) # get Qs
    err_props = map(qual2prop, quals) # convert Qs to propabilities
    avg_err_prop = sum(err_props) / len(qual_str) # calculate average propability
    return prop2qual(avg_err_prop)
# end def get_read_avg_qual


def fastq_packets(fastq, packet_size, num_done_seqs, max_seq_len=None):

    how_to_open = OPEN_FUNCS[ is_gzipped(fastq) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fastq) ]

    with how_to_open(fastq) as fastq_file:

        for _ in range(int(num_done_seqs * FASTQ_LINES_PER_READ)):
            fastq_file.readline()
        # end for

        packet = ""
        qual_dict = dict()
        eof = False

        while not eof:

            for i in range(packet_size):

                read_id = fmt_func(fastq_file.readline())

                if read_id == "":
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

            if not max_seq_len is None:
                packet = prune_seqs(packet.strip(), 'l', max_seq_len)
            # end if

            yield {"fasta": packet.strip(), "qual": qual_dict}

            packet = ""
            qual_dict.clear()
        # end while
    # end with

# end def fastq_packets