
from typing import TextIO

from src.Containers.Fasta import Fasta
from src.WriterSystem.FileWriter import FileWriter


class FastaWriter(FileWriter):

    def _write_single_record(self,
                             sec_record : Fasta,
                             out_file_handle : TextIO):
        out_file_handle.write(f'>{sec_record.header}\n')
        out_file_handle.write(f'{sec_record.seq}\n')
    # end def
# end class

