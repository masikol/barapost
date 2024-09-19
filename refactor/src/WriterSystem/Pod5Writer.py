
import os

from pod5 import Writer

from src.Config.config import OUTPUT_DIR

from src.WriterSystem.FileWriter import FileWriter

from src.Containers.Pod5 import Pod5

# TODO: if this labled file do exist throw "Input path already exists. Refusing to overwrite:"

class Pod5Writer(FileWriter):

    def _write_single_record(self,
                             sec_record : Pod5,
                             out_file_handle : Writer):
        out_file_handle.add_read(sec_record.record.to_read())
    # end def

    def _open_new_outfile(self, outfpath : str) -> Writer:
        return Writer(outfpath)
    # end def

    def _get_out_file_path(self, label : str, index : str) -> str:
        return os.path.join(
            OUTPUT_DIR,
            f'{label}_{index}.{self.ext}'
        )
    # end def
# end class
