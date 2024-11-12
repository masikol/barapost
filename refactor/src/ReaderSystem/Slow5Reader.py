
from pyslow5 import Open

from src.ReaderSystem.FileReader import FileReader
from src.Containers.Slow5 import Slow5


class Slow5Reader(FileReader):

    def _check_file_end(self, record : Slow5) -> bool:
        return False
    # end def

    def _read_single_record(self) -> Slow5:
        try:
            record_id = next(self.reader_iteraror)
            record = self.reader.get_read(record_id)
        except StopIteration:
            raise
        # end try

        return Slow5(record = record)
    # end def

    def open(self) -> None:
        self.reader = Open(self.file_path, 'r')
        self.reader_iteraror = iter(self.reader.get_read_ids()[0])
    # end def
# end class
