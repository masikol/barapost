
import os
import logging
from typing import Generator

from src.Containers.SeqRecord import SeqRecord

from src.ReaderSystem.FileReader import FileReader
from src.ReaderSystem.FastaReader import FastaReader
from src.ReaderSystem.FastQReader import FastQReader

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


class ReaderWrapper(object):

    def __init__(self, 
                 file_path : str, 
                 packet_size : int = 1,
                 probing_packet_size : int = -1,
                 mode : str = 'seq_count',
                 max_seq_len : int = -1):
        
        if not os.path.exists(file_path):
            raise FileExistsError(
                f'File "{file_path}" not found.'
            )
        # end if

        self._validate_positive_integer(packet_size, 'packet_size')
        self._validate_positive_integer(probing_packet_size, 'probing_packet_size', allow_negative_one = True)
        self._validate_positive_integer(max_seq_len, 'max_seq_len', allow_negative_one = True)

        if mode not in ('seq_count', 'sum_seq_len'):
            logger.warning(f'Invalid {mode} mode. Set mode to "seq_count".')
            mode = 'seq_count'
        # end if

        if mode == 'seq_count' and max_seq_len != -1:
            logger.warning(
                f'"max_seq_len" parametr avalible only in "sum_seq_len" mode. Ignore "max_seq_len" parametr.'
                )
            max_seq_len = -1
        # end if
        
        path_split = file_path.split('.')

        if path_split[-1] == 'fasta' or path_split[-2] == 'fasta':
            self.reader = FastaReader(file_path = file_path,
                                     _gzip_ = path_split[-2] == 'fasta',
                                     packet_size = packet_size,
                                     probing_packet_size = probing_packet_size,
                                     mode = mode, 
                                     max_seq_len = max_seq_len)
        elif path_split[-1] == 'fastq' or path_split[-2] == 'fastq':
            self.reader = FastQReader(file_path = file_path,
                                     _gzip_ = path_split[-2] == 'fastq',
                                     packet_size = packet_size,
                                     probing_packet_size = probing_packet_size,
                                     mode = mode, 
                                     max_seq_len = max_seq_len)
        else:
            raise ValueError(
                f'Invalid file type. Use fasta(.gz), fastq(.gz), pod5, fast5, blow5 ot slow5 instead of {path_split[-1]}.'
            )
        # end if
    # end def

    def _validate_positive_integer(self, value: int, name: str, allow_negative_one: bool = False):
        if value < 0 and not (allow_negative_one and value == -1):
            raise ValueError(f'"{name}" must be a positive integer. Received "{name}" = {value}.')
        # end if
    # end def

    def __next__(self) -> Generator[list[SeqRecord], None, None]:   
        return self.reader.__next__()
    # end def

    def __enter__(self) -> FileReader:
        self.reader.open()
        return self.reader
    # end def

    def __exit__(self, type, value, traceback):
        self.reader.close()
    # end def
# end class
