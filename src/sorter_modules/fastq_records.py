# -*- coding: utf-8 -*-

from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped

def fastq_records(fq_path):
    """
    :param read_file: file instance of FASTQ file to retrieve sequences from;
    :type fasta_file: _io.TextIOWrapper or gzip.GzipFile;
    :param fmt_func: function from 'FORMATTING_FUNCS' tuple;

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
    # count number of reads in file (4 lines per read)
    num_reads = int(sum(1 for line in how_to_open(fq_path) if line != '') / 4)

    with how_to_open(fq_path) as fq_file:

        for _ in range(num_reads):
            # Read 4 lines of fastq-record
            yield {
                "seq_id": fmt_func(fq_file.readline()),
                "seq": fmt_func(fq_file.readline()),
                "opt_id": fmt_func(fq_file.readline()),
                "qual_line": fmt_func(fq_file.readline())
            }
        # end for
    # end with
# end def read_fastq_record
