
from gzip import open as gzip_open
from typing import Generator

from src.Containers.SeqRecord import SeqRecord

class FileReader:

    def __init__(self, 
                 file_path : str,
                 _gzip_ : bool = False,
                 packet_size : int = 1, 
                 probing_packet_size : int = -1, 
                 mode : str = 'seq_count',
                 max_seq_len : int = -1):
        
        self.file_path = file_path
        self.packet_size = packet_size
        self.probing_packet_size = probing_packet_size
        self.mode = mode
        self.max_seq_len = max_seq_len

        sum_mode_map = {
            -1: lambda seq : len(seq),
        }
        default_sum_func =  lambda seq : min(self.max_seq_len, len(seq))

        self._open_func = gzip_open if _gzip_ else open
        self._sum_func = sum_mode_map.get(max_seq_len, default_sum_func)
        self._total_records = 0
    # end def


    def _read_single_record(self) -> SeqRecord:
        raise NotImplementedError()
    # end def


    def _check_file_end(self, file : SeqRecord) -> bool:
        raise NotImplementedError()
    # end def

    def _common_generator(self) -> list[SeqRecord]:
        packet = []
        current_sum = 0

        common_condition = lambda condition: condition() and (
            self.probing_packet_size == -1 or self._total_records < self.probing_packet_size
            )

        seq_count_condition = lambda: len(packet) < self.packet_size
        sum_seq_len_condition = lambda: current_sum < self.packet_size
        
        seq_count_condition = lambda: len(packet) < self.packet_size
        sum_seq_len_condition = lambda: current_sum < self.packet_size

        condition = seq_count_condition if self.mode == 'seq_count' else sum_seq_len_condition

        final_condition = lambda: common_condition(condition)

        while final_condition():
            file = self._read_single_record()
            if self._check_file_end(file):
                if packet:
                    return packet
                # end if
                raise StopIteration
            # end if
            packet.append(file)

            if self.mode == 'sum_seq_len':
                current_sum += self._sum_func(file.seq)
            # end if
            
            self._total_records += 1
        # end while

        return packet
    # end def


    def __next__(self) -> Generator[list[SeqRecord], None, None]:   
        if not self.reader:
            raise FileNotFoundError
        # end if

        if self._total_records >= self.probing_packet_size and self.probing_packet_size != -1:
            raise StopIteration
        # end if

        return self._common_generator()
    # end def


    def open(self) -> None:
        self.reader = self._open_func(self.file_path, mode = 'rt')
    # end def


    def close(self) -> None:
        if self.reader:
            self.reader.close()
        # end if
    # end def 

