
import os


FASTA_EXTENSIONS = {
    'fasta',
    'fa',
    'fna',
    'fsa',
    'fasta_nt',
    'fa_nt',
    'fna_nt',
    'fsa_nt',
}

FASTQ_EXTENSIONS = {
    'fastq',
    'fq'
}


def is_fasta(file_path : str) -> bool:
    if is_gzipped(file_path):
        file_path = file_path[:-3]
    # end if
    extension = get_file_extension(file_path)
    return extension.lower() in FASTA_EXTENSIONS
# end def


def is_fastq(file_path : str) -> bool:
    if is_gzipped(file_path):
        file_path = file_path[:-3]
    # end if
    extension = get_file_extension(file_path)
    return extension.lower() in FASTQ_EXTENSIONS
# end def


def is_gzipped(file_path : str) -> bool:
    # TODO: check GZip file validity
    return file_path.endswith('.gz')
# end def


def get_file_extension(file_path : str) -> str:
    return file_path.split('.')[-1]
# end def


def get_hts_file_type(file_path : str) -> str:
    # hts: High Throughput Sequencing
    if is_gzipped(file_path):
        file_path = file_path[:-3]
    # end if
    return get_file_extension(file_path)
# end def


def is_fast5(file_path : str) -> bool:
    extension = get_file_extension(file_path)
    return extension.lower() == 'fast5'
# end def


def is_pod5(file_path : str) -> bool:
    extension = get_file_extension(file_path)
    return extension.lower() == 'pod5'
# end def


def is_blow5(file_path : str) -> bool:
    extension = get_file_extension(file_path)
    return extension.lower() == 'blow5'
# end def


def is_slow5(file_path : str) -> bool:
    extension = get_file_extension(file_path)
    return extension.lower() == 'slow5'
# end def
