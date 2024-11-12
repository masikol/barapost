
from pyslow5 import Open

from src.ReaderSystem.FileReader import FileReader
from src.Containers.Blow5 import Blow5


class Blow5Reader(FileReader):

    def _check_file_end(self, record : Blow5) -> bool:
        return False
    # end def

    def _read_single_record(self) -> Blow5:
        try:
            record_id = next(self.reader_iteraror)
            record = self.reader.get_read(record_id)
        except StopIteration:
            raise
        # end try

        return Blow5(record = record)
    # end def

    def open(self) -> None:
        self.reader = Open(self.file_path, 'r')
        self.reader_iteraror = iter(self.reader.get_read_ids()[0])
    # end def
# end class
