
import gzip
from io import TextIOWrapper
from abc import ABC, abstractmethod
from typing import Generator, Callable, Sequence, MutableSequence

import src.filesystem as fs
from src.Containers.SeqRecord import SeqRecord


class FileReader(ABC):

    def __init__(self,
                 file_paths : Sequence[str],
                 packet_size : int = 1,
                 probing_batch_size : int = -1,
                 mode : str = 'seq_count',
                 max_seq_len : int = -1):

        self.file_paths = file_paths
        # TODO: catch StopIteration
        # Or we will handle this at the arg parsing stage?
        self.curr_file_path = next(iter(file_paths))
        self.curr_file_i = 0

        self.packet_size = packet_size
        self.probing_batch_size = probing_batch_size
        self.mode = mode
        self.max_seq_len = max_seq_len

        self._packet = []
        self._sum_seq_len_read = 0
        self._n_records_read_total = 0

        if self.mode == 'seq_count':
            self._make_packet = self._make_seq_count_packet
        elif self.mode == 'sum_seq_len':
            self._make_packet = self._make_sum_seq_len_packet
        else:
            raise ValueError(
                f'Indalid mode: `{self.mode}`. Allowed nodes: `seq_count`, `sum_seq_len`.'
            )
        # end if

        if max_seq_len == -1:
            self._increment_sum = lambda seq : len(seq)
        elif max_seq_len > 0:
            self._increment_sum = lambda seq : min(self.max_seq_len, len(seq))
        else:
            # TODO: this won't work is str (or sth, not int/float) is passed as max_seq_len
            # Anyway, this validation will be done at the arg parsing stage
            raise ValueError(
                f'Indalid max_seq_len: `{max_seq_len}`. It must be a positive integer or -1.'
            )
        # end if
    # end def

    @abstractmethod
    def _read_single_record(self) -> SeqRecord:
        raise NotImplementedError()
    # end def

    @abstractmethod
    def _check_file_end(self, record : SeqRecord) -> bool:
        raise NotImplementedError()
    # end def

    def _make_conditional_packet(self,
                                 condition : Callable) -> MutableSequence[SeqRecord]:

        while condition():
            try:
                record = self._read_single_record()
            except StopIteration:
                self._switch_to_next_input_file()
                record = self._read_single_record()
            # end try

            if self._check_file_end(record):
                if len(self._packet) != 0:
                    packet = self._packet
                    self._reset_packet()
                    return packet
                # end if
                self._switch_to_next_input_file()
                record = self._read_single_record()
            # end if
            self._packet.append(record)

            if self.mode == 'sum_seq_len':
                self._sum_seq_len_read += self._increment_sum(record.seq)
            # end if

            self._n_records_read_total += 1
        # end while

        packet = self._packet
        self._reset_packet()

        return packet
    # end def

    def _switch_to_next_input_file(self):
        # TODO: this won't work if self.file_paths is a generator
        #   It won't be a generator, anyway, so let it be so
        self.curr_file_i += 1
        if self.curr_file_i >= len(self.file_paths):
            raise StopIteration
        # end if
        self.curr_file_path = self.file_paths[self.curr_file_i]
        self.reader = self._open_gzipwise(self.curr_file_path)
    # end def

    def _reset_packet(self):
        self._packet = []
        self._sum_seq_len_read = 0
    # end def


    def _make_seq_count_packet(self) -> MutableSequence[SeqRecord]:
        return self._make_conditional_packet(
            condition=self._seq_count_stop_condition
        )
    # end def

    def _make_sum_seq_len_packet(self) -> MutableSequence[SeqRecord]:
        return self._make_conditional_packet(
            condition=self._sum_seq_len_stop_condition
        )
    # end def

    def _seq_count_stop_condition(self) -> bool:
        return len(self._packet) < self.packet_size \
               and self._common_stop_condition()
    # end def

    def _sum_seq_len_stop_condition(self) -> bool:
        return self._sum_seq_len_read < self.packet_size \
               and self._common_stop_condition()
    # end def

    def _common_stop_condition(self) -> bool:
        return self.probing_batch_size == -1 \
               or self._n_records_read_total < self.probing_batch_size
    # end def


    # TODO: add return type hint
    def __iter__(self):
        return self
    # end def

    def __next__(self) -> Generator[MutableSequence[SeqRecord], None, None]:
        stop = self.probing_batch_size != -1 \
               and self._n_records_read_total >= self.probing_batch_size
        if stop:
            raise StopIteration
        # end if

        return self._make_packet()
    # end def

    def open(self) -> None:
        self.reader = self._open_gzipwise(self.curr_file_path)
    # end def

    def _open_gzipwise(self, infile_path : str) -> TextIOWrapper:
        if fs.is_gzipped(infile_path):
            return gzip.open(infile_path, mode='rt')
        else:
            return open(infile_path, mode='rt')
        # end if
    # end def

    def close(self) -> None:
        if self.reader:
            self.reader.close()
        # end if
    # end def
# end class
