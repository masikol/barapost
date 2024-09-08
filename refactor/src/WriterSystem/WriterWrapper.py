
import logging

from src.WriterSystem.FileWriter import FileWriter
from src.WriterSystem.FastaWriter import FastaWriter
from src.WriterSystem.FastQWriter import FastQWriter
from src.WriterSystem.Pod5Writer import Pod5Writer
from src.WriterSystem.Fast5Writer import Fast5Writer

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


class WriterWrapper(object):
    def __init__(self, _gzip_ : bool, n_max_out : int, _type_ : str):
        if _type_ == 'fasta':
            self.writer = FastaWriter(_gzip_, n_max_out, 'fasta')
        elif _type_ == 'fastq':
            self.writer = FastQWriter(_gzip_, n_max_out, 'fastq')
        elif _type_ == 'pod5':
            self.writer = Pod5Writer(False, n_max_out, 'pod5')
        elif _type_ == 'fast5':
            self.writer = Fast5Writer(False, n_max_out, 'fast5')
        else:
            raise ValueError(
                f'Invalid file type. Use fasta, fastq or pod5 instead of {_type_}.'
            )
        # end if
        
        if _gzip_ and _type_ in ('pod5', 'fast5'):
            logger.warning(f'Cannot gzip {_type_} format!')
        # end if
    # end def

    def __enter__(self) -> FileWriter:
        return self.writer
    # end def

    def __exit__(self, type, value, traceback):
        self.writer.close()
    # end def
# end class
