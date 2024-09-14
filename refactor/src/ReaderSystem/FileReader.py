
from gzip import open as gipz_open
from typing import Generator

from src.Containers.Fasta import Fasta

class FileReader:
    def __init__(self, 
                 file_path : str,
                 _gzip_ : bool = False,
                 packet_size : int = 1, 
                 probing_packet_size : int = -1, 
                 _mode_ = 'seq_count'):
        self.file_path = file_path
        self.packet_size = packet_size
        self.probing_packet_size = probing_packet_size
        self._mode_ = _mode_

        self._open_func = gipz_open if _gzip_ else open
        self._total_records = 0

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
        if not self.reader:
            raise FileNotFoundError
        # end if

        if self._total_records >= self.probing_packet_size or self.probing_packet_size == -1:
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
            self._total_records += 1
        return packet
    # end def

    def open(self) -> None:
        self.reader = self._open_func(self.file_path, mode = 'rt')
    # end def

    def close(self) -> None:
        if self.reader:
            self.reader.close()
        # end if
    # end def 

