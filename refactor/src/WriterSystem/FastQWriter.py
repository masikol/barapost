import os
import glob
import gzip

from ..Containers.FastQContainer import FastQContainer

from ..Config.config import OUTPUT_DIR

class FastQWriter():
    def __init__(self, _gzip_: bool, N_MAX_OUT: int):
        self._gzip_ = _gzip_
        self.N_MAX_OUT = N_MAX_OUT
        self.file_record_count = {}
    #end def

    def get_last_file_index(self, label):
        extension = "fastq.gz" if self._gzip_ else "fastq"
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

    def write(self, seqs : list[FastQContainer]):
        for fastq in seqs:
            self.write_fastq(fastq)
        # end for
    # end def

    def write_fastq(self, fastq : FastQContainer):
        label = fastq.label
        extension = "fastq.gz" if self._gzip_ else "fastq"

        if label not in self.file_record_count:
            last_index = self.get_last_file_index(label)
            last_filename = os.path.join(OUTPUT_DIR, f'{label}_{last_index}.{extension}')
            record_count = self.get_record_count(last_filename)
            self.file_record_count[label] = (last_index, record_count)
        else:
            last_index, record_count = self.file_record_count[label]
            last_filename = os.path.join(OUTPUT_DIR, f'{label}_{last_index}.{extension}')
        #end if

        if record_count >= self.N_MAX_OUT:
            last_index += 1
            last_filename = os.path.join(OUTPUT_DIR, f'{label}_{last_index}.{extension}')
            record_count = 0
        #end if

        open_func = gzip.open if self._gzip_ else open

        with open_func(last_filename, 'at') as f:
            f.write(f"@{fastq.file.header}\n")
            f.write(f"{fastq.file.seq}\n")
            f.write(f"{fastq.file.plus_line}\n")
            f.write(f"{fastq.file.quality}\n")
        # end with

        record_count += 1
        self.file_record_count[label] = (last_index, record_count)
    # end def

# end class


