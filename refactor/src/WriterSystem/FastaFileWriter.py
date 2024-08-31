import os
import glob
import gzip

from ..Containers.Fasta import Fasta
# from ..Containers.FastQ import FastQ

from ..Config.config import OUTPUT_DIR

class FastaFileWriter():
    def __init__(self, _gzip_: bool, n_max_out: int):
        self._gzip_ = _gzip_
        self.n_max_out = n_max_out
        self.file_record_count = {}
    # end def

    def get_last_file_index(self, label):
        extension = "fasta.gz" if self._gzip_ else "fasta"
        pattern = os.path.join(OUTPUT_DIR, f"{label}_*.{extension}")
        files = glob.glob(pattern)
        
        if not files:
            return 0
        # end if

        indices = []
        for f in files:
            basename = os.path.basename(f)
            try:
                index_str = basename.split('_')[-1].split('.')[0]
                index = int(index_str)
                indices.append(index)
            except ValueError:
                pass
            # end try

        return max(indices) if indices else 0
    # end def 

    def get_record_count(self, label):
        if label not in self.file_record_count:
            return 0
        # end if
        return self.file_record_count[label]
    # end def

    def write(self, seqs: list):
        for fasta in seqs:
            self.write_fasta(fasta)
        # end for
    # end def 

    def write_fasta(self, fasta: Fasta):
        label = fasta.label
        extension = "fasta.gz" if self._gzip_ else "fasta"

        if label not in self.file_record_count:
            last_index = self.get_last_file_index(label)
            
            last_filename = os.path.join(OUTPUT_DIR, f'{label}_{last_index}.{extension}')
            record_count = self.get_record_count(label)
            self.file_record_count[label] = (last_index, record_count)
        else:
            last_index, record_count = self.file_record_count[label]
            last_filename = os.path.join(OUTPUT_DIR, f'{label}_{last_index}.{extension}')
        # end if

        if record_count >= self.n_max_out:
            last_index += 1
            last_filename = os.path.join(OUTPUT_DIR, f'{label}_{last_index}.{extension}')
            record_count = 0
        # end if

        open_func = gzip.open if self._gzip_ else open

        with open_func(last_filename, 'at') as f:
            f.write(f">{fasta.fastaHeader}\n")
            f.write(f"{fasta.seq}\n")
        # end with

        record_count += 1
        self.file_record_count[label] = (last_index, record_count)
    # end def
# end class

