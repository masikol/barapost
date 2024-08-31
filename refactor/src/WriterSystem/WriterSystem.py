from ..WriterSystem.FastaFileWriter import FastaFileWriter
from ..WriterSystem.FastQWriter import FastQWriter

class WriterSystem():
    def __init__(self, _gzip_ : bool, n_max_out : int, _type_ : str):
        self._gzip_ = _gzip_
        self.n_max_out = n_max_out
        if _type_ not in ['fasta', 'fastq']:
            raise ValueError(f'Invalid file type! Use fasta or fastq instead of {_type_}!')
        # end if
        self._type_ = _type_

        if _type_ == 'fasta':
            self.writer = FastaFileWriter(_gzip_, n_max_out)
        elif _type_ == 'fastq':
            self.writer = FastQWriter(_gzip_, n_max_out)
    # end def

    def write(self, seqs : list):
        self.writer.write(seqs)
    # end def
#end class

