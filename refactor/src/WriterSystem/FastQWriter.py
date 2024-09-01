
from typing import TextIO

from src.Containers.FastQ import FastQ
from src.WriterSystem.FileWriter import FileWriter


class FastQWriter(FileWriter):

    def _write_single_record(self,
                             sec_record : FastQ,
                             out_file_handle : TextIO):
        out_file_handle.write(f'@{sec_record.header}\n')
        out_file_handle.write(f'{sec_record.seq}\n')
        out_file_handle.write(f'{sec_record.plus_line}\n')
        out_file_handle.write(f'{sec_record.quality}\n')
    # end def
# end class
