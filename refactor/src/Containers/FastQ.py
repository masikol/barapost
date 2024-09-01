
import math
import statistics


class FastQ:
    def __init__(self,
                 header : str,
                 seq : str,
                 plus_line : str,
                 quality : str,
                 offset : int = 33):
        self.header = header
        self.seq = seq
        self.plus_line = plus_line
        self.quality = quality
        if offset not in (33, 64):
            # TODO: print info messages with logging
            print(f'Unexpected offset = {offset}! Offset was set to 33.')
            self.offset = 33
        else:
            self.offset = offset
        # end if
    # end def

    def average_quality(self) -> float:
        avg_error_prob = statistics.mean(
            map(self._phred_char_to_pe, self.quality)
        )
        return self._pe_to_Q(avg_error_prob)
    # end def

    def _phred_char_to_pe(self, char : str) -> float:
        # pe is error probability
        Q = ord(char) - self.offset
        return 10 ** (-Q / 10)
    # end def

    def _pe_to_Q(self, error_prob : float):
        return -10 * math.log10(error_prob)
    # end def
# end class
