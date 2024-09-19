
import gzip
from abc import ABC, abstractmethod
from typing import Generator, Callable, MutableSequence

from src.Containers.SeqRecord import SeqRecord

CURRENT_SUM_IDX = 0


class FileReader(ABC):

    def __init__(self,
                 file_path : str,
                 _gzip_ : bool = False,
                 packet_size : int = 1,
                 probing_batch_size : int = -1,
                 mode : str = 'seq_count',
                 max_seq_len : int = -1):

        self.file_path = file_path
        self.packet_size = packet_size
        self.probing_batch_size = probing_batch_size
        self.mode = mode
        self.max_seq_len = max_seq_len

        self._total_records = 0
        self._open_func = gzip.open if _gzip_ else open

        generators_map = {
            'seq_count' : self._seq_count_generator,
            'sum_seq_len' : self._sum_seq_len_generator
        }
        self.generator_func = generators_map[mode]

        if max_seq_len == -1:
            self._sum_func = lambda seq : len(seq)
        else:
            self._sum_func = lambda seq : min(self.max_seq_len, len(seq))
        # end if

        self.common_condition = lambda condition: condition() and (
            self.probing_batch_size == -1 or self._total_records < self.probing_batch_size
        )
    # end def

    @abstractmethod
    def _read_single_record(self) -> SeqRecord:
        raise NotImplementedError()
    # end def

    @abstractmethod
    def _check_file_end(self, record : SeqRecord) -> bool:
        raise NotImplementedError()
    # end def

    def _common_generator(self,
                          condition : Callable,
                          packet : MutableSequence,
                          current_sum : MutableSequence[int]) -> MutableSequence[SeqRecord]:

        while condition():
            record = self._read_single_record()
            if self._check_file_end(record):
                if packet:
                    return packet
                # end if
                raise StopIteration
            # end if
            packet.append(record)

            if self.mode == 'sum_seq_len':
                current_sum[CURRENT_SUM_IDX] += self._sum_func(record.seq)
            # end if

            self._total_records += 1
        # end while

        return packet
    # end def

    def _seq_count_generator(self,
                             packet : MutableSequence,
                             current_sum : MutableSequence[int]) -> MutableSequence[SeqRecord]:
        seq_count_condition = lambda: len(packet) < self.packet_size
        final_condition = lambda: self.common_condition(seq_count_condition)
        return self._common_generator(
            condition = final_condition,
            packet = packet,
            current_sum = current_sum
        )
    # end def

    def _sum_seq_len_generator(self,
                               packet : MutableSequence,
                               current_sum : MutableSequence[int]) -> MutableSequence[SeqRecord]:
        sum_seq_len_condition = lambda: current_sum[CURRENT_SUM_IDX] < self.packet_size
        final_condition = lambda: self.common_condition(sum_seq_len_condition)
        return self._common_generator(
            condition = final_condition,
            packet = packet,
            current_sum = current_sum
        )
    # end def

    def __next__(self) -> Generator[MutableSequence[SeqRecord], None, None]:
        if self._total_records >= self.probing_batch_size and self.probing_batch_size != -1:
            raise StopIteration
        # end if

        packet = []
        current_sum = [0] # Cause int immutable

        return self.generator_func(packet = packet, current_sum = current_sum)
    # end def

    def open(self) -> None:
        self.reader = self._open_func(self.file_path, mode = 'rt')
    # end def

    def close(self) -> None:
        if self.reader:
            self.reader.close()
        # end if
    # end def
# end class
