
import logging

from src.ReaderSystem.FileReader import FileReader

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


class ReaderWrapper(object):
    def __init__(self, 
                 file_path : str, 
                 packet_size : int = 1,
                 probing_packet_size : int = 1,
                 _mode_ : str = 'seq_count'):
        
        path_split = file_path.split('.')

        if path_split[-1] == 'fasta' or path_split[-2] == 'fasta':
            self.reader = FileReader(file_path = file_path,
                                     _gzip_ = path_split[-2] == 'fasta',
                                     packet_size = packet_size,
                                     probing_packet_size = probing_packet_size,
                                     _mode_ = _mode_)
        else:
            raise ValueError(
                f'Invalid file type. Use fasta, fastq, pod5, fast5, blow5 ot slow5 instead of {file_path}.'
            )
        # end if
    # end def

    def __next__(self):   
        return self.reader.next()
    # end def

    def __enter__(self) -> FileReader:
        self.reader.open()
        return self.reader
    # end def

    def __exit__(self, type, value, traceback):
        self.reader.close()
    # end def
# end class
