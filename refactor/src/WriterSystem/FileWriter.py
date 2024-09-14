
import os
import gzip
from typing import TextIO, MutableSequence
from abc import ABC, abstractmethod

from src.Config.config import OUTPUT_DIR

from src.Containers.ClassifContainer import ClassifContainer

OUT_FILE_HANDLE_IDX = 0
# CURR_INDEX_IDX      = 1 # not used for now but we leave it here
RECORD_COUNT_IDX    = 2


class FileWriter(ABC):

    @abstractmethod
    def _write_single_record(self, sec_record, outfpath):
        raise NotImplementedError()
    # end def


    def __init__(self, _gzip_: bool, n_max_out: int, ext = 'fasta'):
        self._gzip_ = _gzip_
        self.open_func = gzip.open if self._gzip_ else open
        self.n_max_out = n_max_out
        self.ext = ext
        self.record_chunk_counter = dict()
    # end def

    def write(self, classif_records : MutableSequence[ClassifContainer]):
        for cr in classif_records:
            out_file_handle = self._get_out_file_handle(cr.label)
            self._write_single_record(cr.record, out_file_handle)
            self._update_record_count(cr.label)
        # end for
    # end def

    def _get_out_file_handle(self, label) -> str:
        if label in self.record_chunk_counter:
            out_file_handle = self.record_chunk_counter[label][
                OUT_FILE_HANDLE_IDX
            ]
        else:
            curr_index, record_count = 0, 0
            outfpath = self._get_out_file_path(label, curr_index)
            out_file_handle = self._open_new_outfile(outfpath)
            self.record_chunk_counter[label] = [
                out_file_handle,
                curr_index,
                record_count
            ]
        # end try
        return out_file_handle
    # end def

    def _open_new_outfile(self, outfpath : str) -> TextIO:
        return self.open_func(outfpath, 'wt')
    # end def

    def _update_record_count(self, label):
        outfile_handle, curr_index, record_count = self.record_chunk_counter[label]
        record_count += 1
        self.record_chunk_counter[label][RECORD_COUNT_IDX] = record_count

        if record_count >= self.n_max_out:
            outfile_handle.close()
            curr_index += 1
            record_count = 0
            outfpath = self._get_out_file_path(label, curr_index)
            outfile_handle = self._open_new_outfile(outfpath)

            self.record_chunk_counter[label] = [
                outfile_handle,
                curr_index,
                record_count
            ]
        # end if
    # end def

    def _get_out_file_path(self, label : str, index : str) -> str:
        extension = self.ext
        if self._gzip_:
            extension = f'{extension}.gz'
        # end if
        return os.path.join(
            OUTPUT_DIR,
            f'{label}_{index}.{extension}'
        )
    # end def

    def close(self):
        for v_tuple in self.record_chunk_counter.values():
            v_tuple[OUT_FILE_HANDLE_IDX].close()
        # end for
    # end def
# end class
