
from h5py import File

from src.Containers.Fast5 import Fast5
from src.simplify_read_id import simplify_read_id
from src.ReaderSystem.FileReader import FileReader


class Fast5Reader(FileReader):

    def _check_file_end(self, record : Fast5) -> bool:
        return False
    # end def


    def _read_single_record(self) -> Fast5:
        try:
            read_id = next(self.reader_iterator)
        except StopIteration:
            raise
        # end try
        return Fast5(file_handle = self.reader, read_id = read_id)
    # end def

    def _detect_single_fast5(self):
        return 'Raw' in self.reader.keys()
    # end def


    def open(self) -> None:
        self.reader = File(self.file_path, 'r')
        if self._detect_single_fast5():
            read_id = 'read_{}'.format(
                simplify_read_id(self.reader.filename) # read id is in the filename
            )
            self.reader_iterator = iter( (read_id,) )
        else:
            self.reader_iterator = iter(self.reader)
        # end if
    # end def


    # TODO: Manualy close it after writing!!!
    def close(self) -> None: # Do not close FAST5 files till write it
        pass
    # end def
# end class
