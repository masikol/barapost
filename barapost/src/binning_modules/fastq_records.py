# -*- coding: utf-8 -*-
# Module defines generator that yields records retrieved from fastq files.

from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped

def fastq_records(fq_path):
    """
    Generator yields records retrieved from fasta files.

    :param read_file: file instance of FASTQ file to retrieve sequences from;
    :type fasta_file: _io.TextIOWrapper or gzip.GzipFile;

    Returns dictionary of the following structure:
    {
        "seq_id": ID_of_sequence,
        "seq": sequence_itself,
        "opt_id": the_third_line,
        "qual_line": quality_line
    }
    """

    how_to_open = OPEN_FUNCS[ is_gzipped(fq_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_path) ]

    with how_to_open(fq_path) as fq_file:

        eof = False

        while not eof:

            seq_id = fmt_func(fq_file.readline())
            if seq_id != "":
                yield {
                    "seq_id": seq_id,
                    "seq": fmt_func(fq_file.readline()),
                    "opt_id": fmt_func(fq_file.readline()),
                    "qual_line": fmt_func(fq_file.readline())
                }
            else:
                eof = True
            # end if
        # end while
    # end with
# end def read_fastq_record
