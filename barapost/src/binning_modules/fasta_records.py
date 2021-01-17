# -*- coding: utf-8 -*-
# Module defines generator that yields records retrieved from fasta files.

from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped

def fasta_records(fa_path):
    # Generator yields records retrieved from fasta files.
    #
    # :param fasta_file: file instance of FASTA file to retrieve sequences from;
    # :type fasta_file: _io.TextIOWrapper or gzip.GzipFile;
    # Returns dictionary of the following structure:
    # {
    #     "seq_id": ID_of_sequence,
    #     "seq": sequence_itself
    # }

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
# end def read_fasta_record
