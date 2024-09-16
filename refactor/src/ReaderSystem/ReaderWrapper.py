
import os
import logging
from typing import Generator

from src.Containers.SeqRecord import SeqRecord

from src.ReaderSystem.FileReader import FileReader
from src.ReaderSystem.FastaReader import FastaReader
from src.ReaderSystem.FastQReader import FastQReader
from src.ReaderSystem.Fast5Reader import Fast5Reader
from src.ReaderSystem.Pod5Reader import Pod5Reader
from src.ReaderSystem.Slow5Reader import Slow5Reader

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

        self.file_path = file_path
        self.packet_size = packet_size
        self.probing_packet_size = probing_packet_size
        self.mode = mode
        self.max_seq_len = max_seq_len

        path_split = self.file_path.split('.')

        self._validate_positive_integer(self.packet_size, 'packet_size')
        self._validate_positive_integer(self.probing_packet_size, 'probing_packet_size', allow_negative_one = True)
        self._validate_positive_integer(self.max_seq_len, 'max_seq_len', allow_negative_one = True)

        if self.mode == 'seq_count' and self.max_seq_len != -1:
            logger.warning(
                f'"max_seq_len" parametr avalible only in "sum_seq_len" mode. Ignore "max_seq_len" parametr.'
                )
            self.max_seq_len = -1
        # end if

        if path_split[-1] in ('fast5', 'pod5', 'slow5'):
            for name, value, allowed_value in zip(
                    ('packet_size', 'probing_packet_size', 'mode', 'max_seq_len'),
                    (packet_size, probing_packet_size, mode, max_seq_len),
                    (1, -1, 'seq_count', -1)):
                self._validate_input_parametrs(path_split[-1], value, allowed_value, name)
            # end for
        # end if

        if self.mode not in ('seq_count', 'sum_seq_len'):
            logger.warning(f'Invalid {self.mode} mode. Set mode to "seq_count".')
            self.mode = 'seq_count'
        # end if

        if path_split[-1] == 'fasta' or path_split[-2] == 'fasta':
            self.reader = FastaReader(file_path = self.file_path,
                                     _gzip_ = path_split[-2] == 'fasta',
                                     packet_size = self.packet_size,
                                     probing_packet_size = self.probing_packet_size,
                                     mode = self.mode, 
                                     max_seq_len = self.max_seq_len)
        elif path_split[-1] == 'fastq' or path_split[-2] == 'fastq':
            self.reader = FastQReader(file_path = self.file_path,
                                     _gzip_ = path_split[-2] == 'fastq',
                                     packet_size = self.packet_size,
                                     probing_packet_size = self.probing_packet_size,
                                     mode = self.mode, 
                                     max_seq_len = self.max_seq_len)
        elif path_split[-1] == 'fast5':
            self.reader = Fast5Reader(file_path = self.file_path,
                                     _gzip_ = False,
                                     packet_size = self.packet_size,
                                     probing_packet_size = self.probing_packet_size,
                                     mode = self.mode, 
                                     max_seq_len = self.max_seq_len)
        elif path_split[-1] == 'pod5':
            self.reader = Pod5Reader(file_path = self.file_path,
                                     _gzip_ = False,
                                     packet_size = self.packet_size,
                                     probing_packet_size = self.probing_packet_size,
                                     mode = self.mode, 
                                     max_seq_len = self.max_seq_len)
        elif path_split[-1] == 'slow5':
            self.reader = Slow5Reader(file_path = self.file_path,
                                     _gzip_ = False,
                                     packet_size = self.packet_size,
                                     probing_packet_size = self.probing_packet_size,
                                     mode = self.mode, 
                                     max_seq_len = self.max_seq_len)
        else:
            raise ValueError(
                f'Invalid file type. Use fasta(.gz), fastq(.gz), pod5, fast5, blow5 or slow5 instead of {path_split[-1]}.'
            )
        # end if
    # end def

    def _validate_positive_integer(self, value : int, name : str, allow_negative_one : bool = False):
        if value < 0 and not (allow_negative_one and value == -1):
            raise ValueError(f'"{name}" must be a positive integer. Received "{name}"={value}.')
        # end if
    # end def

    def _validate_input_parametrs(self, 
                                  file_type : str, 
                                  value : int,
                                  allowed_value : int,
                                  name : str) -> bool:
        if value != allowed_value:
            logger.warning(
                f'"{name}" parameter is not supported with {file_type} format. Setting "{name}" to "{allowed_value}".'
            )
            setattr(self, name, allowed_value)
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
