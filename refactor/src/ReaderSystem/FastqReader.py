
from src.ReaderSystem.FileReader import FileReader
from src.Containers.Fastq import Fastq


class FastqReader(FileReader):

    def _check_file_end(self, record : Fastq) -> bool:
        return record.header == ''
    # end def

    def _read_single_record(self) -> Fastq:
        header    = self.reader.readline().strip('\r\n')
        seq       = self.reader.readline().strip('\r\n')
        plus_line = self.reader.readline().strip('\r\n')
        quality   = self.reader.readline().strip('\r\n')

        return Fastq(
            header=header,
            seq=seq,
            plus_line=plus_line,
            quality=quality
        )
    # end def
# end class
