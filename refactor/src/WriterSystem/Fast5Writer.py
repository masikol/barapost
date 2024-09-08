
import os
from h5py import File

from src.Config.config import OUTPUT_DIR

from src.Containers.Fast5 import Fast5

from src.WriterSystem.FileWriter import FileWriter


class Fast5Writer(FileWriter):

    def _write_single_record(self,
                             sec_record : Fast5,
                             out_file_handle : File):
        sec_record.out_file_handle.copy(sec_record.read_uuid, out_file_handle)
    # end def

    def _get_out_file_path(self, label : str, index : str) -> str:
        return os.path.join(
            OUTPUT_DIR,
            f'{label}_{index}.{self.ext}'
        )
    # end def

    def _open_new_outfile(self, outfpath : str) -> File:
        return File(outfpath, 'w')
    # end def

# end class

