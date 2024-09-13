
from src.Containers.Fasta import Fasta

from typing import Generator

class FileReader:
    def __init__(self, 
                 file_path : str,
                 packet_size : int = 1, 
                 probing_packet_size : int = -1, 
                 _mode_ = 'seq_count'):
        self.file_path = file_path
        self.packet_size = packet_size
        self.probing_packet_size = probing_packet_size
        self._mode_ = _mode_
        self.total_records = 0

        # raise NotImplementedError()
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
            if line.startswith(">"):
                self.reader.seek(pos)
                break
            seq_lines.append(line)
            # end if
        # end while

        seq = ''.join(seq_lines)
        return Fasta(header, seq)
    # end def


    def _check_file_end(self, file : Fasta) -> bool:
        return file.header == ''
    # end def


    def __next__(self) -> Generator[list[Fasta], None, None]:   
        return self.next()
    # end def

    def next(self) -> Generator[list[Fasta], None, None]:
        if self.total_records >= self.probing_packet_size or self.probing_packet_size == -1:
            raise StopIteration
        # end if
        packet = []
        for _ in range(self.packet_size):
            file = self._read_single_record()
            if self._check_file_end(file):
                if packet:
                    return packet
                # end if
                raise StopIteration
            packet.append(file)
            self.total_records += 1
        return packet
    # end def

    def open(self) -> None:
        self.reader = open(self.file_path, mode = 'r')
    # end def

    def close(self) -> None:
        self.reader.close()
    # end def 

