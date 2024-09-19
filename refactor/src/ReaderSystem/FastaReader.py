
from src.ReaderSystem.FileReader import FileReader
from src.Containers.Fasta import Fasta


class FastaReader(FileReader):

    def _check_file_end(self, record : Fasta) -> bool:
        return record.header == ''
    # end def

    def _read_single_record(self) -> Fasta:
        header = self.reader.readline().strip('\r\n')
        seq_lines = []

        while True:
            pos = self.reader.tell()
            line = self.reader.readline().strip('\r\n')
            if line == '':
                break
            # end if
            if line.startswith('>'):
                self.reader.seek(pos)
                break
            # end if
            seq_lines.append(line)
        # end while

        seq = ''.join(seq_lines)
        return Fasta(header = header, seq = seq)
    # end def
# end class
