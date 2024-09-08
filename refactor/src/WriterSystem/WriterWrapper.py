
from src.WriterSystem.FileWriter import FileWriter
from src.WriterSystem.FastaWriter import FastaWriter
from src.WriterSystem.FastQWriter import FastQWriter
from src.WriterSystem.Pod5Writer import Pod5Writer


class WriterWrapper(object):
    def __init__(self, _gzip_ : bool, n_max_out : int, _type_ : str):
        if _type_ == 'fasta':
            self.writer = FastaWriter(_gzip_, n_max_out, 'fasta')
        elif _type_ == 'fastq':
            self.writer = FastQWriter(_gzip_, n_max_out, 'fastq')
        elif _type_ == 'pod5':
            self.writer = Pod5Writer(False, n_max_out, 'pod5')
        else:
            raise ValueError(
                f'Invalid file type. Use fasta or fastq instead of {_type_}.'
            )
        # end if
    # end def

    def __enter__(self) -> FileWriter:
        return self.writer
    # end def

    def __exit__(self, type, value, traceback):
        self.writer.close()
    # end def
# end class
