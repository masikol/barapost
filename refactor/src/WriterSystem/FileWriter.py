import os
import gzip

from ..Containers.Fasta import Fasta
from ..Containers.FastQ import FastQ

from ..Config.config import OUTPUT_DIR

class FileWriter:
    def __init__(self, gzip: bool, file_type: str, N_MAX_OUT: int):
        self.gzip = gzip
        self.file_type = file_type
        self.N_MAX_OUT = N_MAX_OUT
        self.file_record_count = {}
    #end def

    def get_last_file_index(self, label):
        index = 0
        while True:
            extension = f"{self.file_type.lower()}.gz" if self.gzip else self.file_type.lower()
            filename = f"{OUTPUT_DIR}/{label}_{index}.{extension}" if index > 0 else f"{OUTPUT_DIR}/{label}.{extension}"
            if not os.path.exists(filename):
                break
            #end if
            index += 1
        #nd while
        return index - 1 if index > 0 else 0
    #end def

    def get_record_count(self, filename):
        if not os.path.exists(filename):
            return 0
        #end if
        with open(filename, 'rt') if not self.gzip else gzip.open(filename, 'rt') as f:
            if self.file_type == 'FASTA':
                return sum(1 for line in f if line.startswith('>'))
            elif self.file_type == 'FASTQ':
                return sum(1 for i, line in enumerate(f) if i % 4 == 0)
            #end if
        #end with
    #end def

    def write(self, seqs: list):
        if self.file_type == 'FASTA':
            for fasta in seqs:
                self.write_fasta(fasta)
            #end for
        elif self.file_type == 'FASTQ':
            for fastq in seqs:
                self.write_fastq(fastq)
            #end for
        #end if
    #end def

    def write_fasta(self, fasta: Fasta):
        label = fasta.label
        if label not in self.file_record_count:
            last_index = self.get_last_file_index(label)
            extension = "fasta.gz" if self.gzip else "fasta"
            last_filename = f"{OUTPUT_DIR}/{label}_{last_index}.{extension}"
            record_count = self.get_record_count(last_filename)
            self.file_record_count[label] = (last_index, record_count)
        else:
            last_index, record_count = self.file_record_count[label]
            extension = "fasta.gz" if self.gzip else "fasta"
            last_filename = f"{OUTPUT_DIR}/{label}_{last_index}.{extension}"
        #end if 

        if record_count >= self.N_MAX_OUT:
            last_index += 1
            last_filename = f"{OUTPUT_DIR}/{label}_{last_index}.{extension}"
            record_count = 0
        #end if

        write_mode = 'ab' if self.gzip else 'a'
        open_func = gzip.open if self.gzip else open

        with open_func(last_filename, write_mode) as f:
            f.write(f">{fasta.fatataHeader}\n".encode())
            f.write(f"{fasta.seq}\n".encode())
        #end with

        record_count += 1
        self.file_record_count[label] = (last_index, record_count)
    #end def

    def write_fastq(self, fastq: FastQ):
        label = fastq.label
        if label not in self.file_record_count:
            last_index = self.get_last_file_index(label)
            extension = "fastq.gz" if self.gzip else "fastq"
            last_filename = f"{OUTPUT_DIR}/{label}_{last_index}.{extension}"
            record_count = self.get_record_count(last_filename)
            self.file_record_count[label] = (last_index, record_count)
        else:
            last_index, record_count = self.file_record_count[label]
            extension = "fastq.gz" if self.gzip else "fastq"
            last_filename = f"{OUTPUT_DIR}/{label}_{last_index}.{extension}"
        #end if

        if record_count >= self.N_MAX_OUT:
            last_index += 1
            last_filename = f"{OUTPUT_DIR}/{label}_{last_index}.{extension}"
            record_count = 0
        #end if

        write_mode = 'ab' if self.gzip else 'a'
        open_func = gzip.open if self.gzip else open

        with open_func(last_filename, write_mode) as f:
            f.write(f"@{fastq.label}\n".encode())
            f.write(f"{fastq.seq}\n".encode())
            f.write(f"{fastq.plus_line}\n".encode())
            f.write(f"{fastq.quality}\n".encode())
        #end with

        record_count += 1
        self.file_record_count[label] = (last_index, record_count)
    #end def

# end class


