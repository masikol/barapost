
import logging

from src.WriterSystem.FileWriter  import FileWriter
from src.WriterSystem.FastaWriter import FastaWriter
from src.WriterSystem.FastqWriter import FastqWriter
from src.WriterSystem.Pod5Writer  import Pod5Writer
from src.WriterSystem.Fast5Writer import Fast5Writer
from src.WriterSystem.Blow5Writer import Blow5Writer
from src.WriterSystem.Slow5Writer import Slow5Writer

# TODO: don't forget to move higher to some config abstraction level
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


class WriterWrapper(object):
    def __init__(self, _gzip_ : bool, n_max_out : int, _type_ : str):

        if _type_.lower() == 'fasta':
            self.writer = FastaWriter(_gzip_, n_max_out, 'fasta')
        elif _type_.lower() == 'fastq':
            self.writer = FastqWriter(_gzip_, n_max_out, 'fastq')
        elif _type_.lower() == 'pod5':
            self.writer = Pod5Writer( False,  n_max_out,  'pod5')
        elif _type_.lower() == 'fast5':
            self.writer = Fast5Writer(False,  n_max_out, 'fast5')
        elif _type_.lower() == 'blow5':
            self.writer = Blow5Writer(False,  n_max_out, 'blow5')
        elif _type_.lower() == 'slow5':
            self.writer = Slow5Writer(False,  n_max_out, 'slow5')
        else:
            raise ValueError(
                f'Invalid `_type_` argument: `{_type_}`. Allowed types: `fasta`, `fastq`, `pod5`, `fast5`, `blow5`, `slow5`.'
            )
        # end if
        
        if _gzip_ and _type_.lower() in ('pod5', 'fast5', 'blow5', 'slow5'):
            logger.warning(f'Cannot gzip {_type_} data. Ignoring output compression.')
        # end if
    # end def

    def __enter__(self) -> FileWriter:
        return self.writer
    # end def

    def __exit__(self, type, value, traceback):
        self.writer.close()
    # end def
# end class
