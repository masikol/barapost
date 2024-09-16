
from h5py import File

from src.ReaderSystem.FileReader import FileReader
from src.Containers.Fast5 import Fast5

class Fast5Reader(FileReader):        

    def _check_file_end(self, record : Fast5) -> bool:
        return False
    # end def

    def _read_single_record(self) -> Fast5:
        try:
            uuid = next(self.reader_iterator)
        except StopIteration:
            raise
        # end try

        return Fast5(out_file_handle = self.reader, read_uuid = uuid)
    # end def

    def open(self) -> None:
        self.reader = File(self.file_path, 'r')
        self.reader_iterator = iter(self.reader)
    # end def

    # TODO: Manualy close it after writing!!!
    def close(self) -> None: # Do not close FAST5 files till write it
        pass

# end class