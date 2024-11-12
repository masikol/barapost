
import math
import logging
import statistics

from src.Containers.SeqRecord import SeqRecord

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)


class Fastq(SeqRecord):

    __slots__ = ('header', 'seq', 'plus_line', 'quality', 'offset')

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
            logger.warning(f'Unexpected offset: `{offset}`. Setting offset to 33.')
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

    def __str__(self):
        seq_concise     = self._get_consice_str(self.seq)
        quality_concise = self._get_consice_str(self.quality)
        return f'''header: {self.header},
seq: {seq_concise},
plus_line: {self.plus_line},
quality: {quality_concise}.\n'''
    # end def

    def __repr__(self):
        seq_concise     = self._get_consice_str(self.seq)
        quality_concise = self._get_consice_str(self.quality)
        return f'''Fastq(
    header={self.header!r}, 
    sequence={seq_concise!r}, 
    plus_line={self.plus_line!r}, 
    quality={quality_concise!r}
)'''
    # end def

    def _get_consice_str(self, string):
        n_chars_show = 30
        n_chars_omitted = len(self.seq) - 2*n_chars_show
        return '{}../{:,}chars/..{}'.format(
            string[:n_chars_show],
            n_chars_omitted,
            string[-n_chars_show:]
        )
    # end def
# end class
