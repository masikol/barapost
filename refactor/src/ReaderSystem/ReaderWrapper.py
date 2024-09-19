
import os
import logging
from typing import Generator, MutableSequence

from src.Containers.SeqRecord import SeqRecord

from src.ReaderSystem.FileReader import FileReader
from src.ReaderSystem.FastaReader import FastaReader
from src.ReaderSystem.FastQReader import FastQReader
from src.ReaderSystem.Fast5Reader import Fast5Reader
from src.ReaderSystem.Pod5Reader import Pod5Reader
from src.ReaderSystem.Slow5Reader import Slow5Reader
from src.ReaderSystem.Blow5Reader import Blow5Reader


# TODO: don't forget to move higher to some config abstraction level
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


class ReaderWrapper(object):

    def __init__(self,
                 file_path : str,
                 packet_size : int = 1,
                 probing_batch_size : int = -1,
                 mode : str = 'seq_count',
                 max_seq_len : int = -1):
        
        if not os.path.isfile(file_path):
            raise FileNotFoundError(
                f'File `{file_path}` does not exist.'
            )
        # end if

        self.file_path          = file_path
        self.packet_size        = packet_size
        self.probing_batch_size = probing_batch_size
        self.mode               = mode
        self.max_seq_len        = max_seq_len

        path_split = self.file_path.split('.')

        self._validate_positive_integer(self.packet_size, 'packet_size')
        self._validate_positive_integer(
            self.probing_batch_size,
            'probing_batch_size',
            allow_minus_one = True
        )
        self._validate_positive_integer(
            self.max_seq_len,
            'max_seq_len',
            allow_minus_one = True
        )

        if self.mode == 'seq_count' and self.max_seq_len != -1:
            logger.warning(
                f'The `max_seq_len` parameter is avalible only in `sum_seq_len` mode.'
            )
            logger.warning('Ignoring the `max_seq_len` parameter.')
            self.max_seq_len = -1
        # end if

        if self.mode not in ('seq_count', 'sum_seq_len'):
            logger.warning(f'Invalid mode: `{self.mode}`. Setting mode to `seq_count`.')
            self.mode = 'seq_count'
        # end if

        # TODO: we should allow multiple variations of common fasta and fastq extentions.
        # fastQ:
        #   fastq, fq
        # fastA:
        #   fasta, fa, fna, fsa, fasta_nt, fa_nt, fna_nt, fsa_nt
        # And all this might be .gz.
        # Variant: make a module `filesystem` and make there functions
        #   `is_fasta`, `is_fastq`, `is_gzipped`
        if path_split[-1] == 'fasta' or path_split[-2] == 'fasta':
            self.reader = FastaReader(
                file_path=self.file_path,
                _gzip_=path_split[-2] == 'fasta',
                packet_size=self.packet_size,
                probing_batch_size=self.probing_batch_size,
                mode=self.mode,
                max_seq_len=self.max_seq_len
            )
        elif path_split[-1] == 'fastq' or path_split[-2] == 'fastq':
            self.reader = FastQReader(
                file_path=self.file_path,
                _gzip_=path_split[-2] == 'fastq',
                packet_size=self.packet_size,
                probing_batch_size=self.probing_batch_size,
                mode=self.mode,
                max_seq_len=self.max_seq_len
            )
        elif path_split[-1] == 'fast5':
            self.reader = Fast5Reader(
                file_path=self.file_path
            )
        elif path_split[-1] == 'pod5':
            self.reader = Pod5Reader(
                file_path=self.file_path
            )
        elif path_split[-1] == 'blow5':
            self.reader = Blow5Reader(
                file_path=self.file_path
            )
        elif path_split[-1] == 'slow5':
            self.reader = Slow5Reader(
                file_path=self.file_path
            )
        else:
            raise ValueError(
                f'Invalid file type: `{path_split[-1]}`. Allowed types: fasta(.gz), fastq(.gz), pod5, fast5, blow5, slow5.'
            )
        # end if
    # end def

    def _validate_positive_integer(self,
                                   value : int,
                                   name : str,
                                   allow_minus_one : bool = False):
        if value < 0 and not (allow_minus_one and value == -1):
            raise ValueError(f'`{name}` must be a positive integer. Received `{name}`=`{value}`')
        # end if
    # end def

    def __next__(self) -> Generator[MutableSequence[SeqRecord], None, None]:
        return next(self.reader)
    # end def

    def __enter__(self) -> FileReader:
        self.reader.open()
        return self.reader
    # end def

    def __exit__(self, type, value, traceback):
        self.reader.close()
    # end def
# end class
