
from src.Containers.SeqRecord import SeqRecord


class Fasta(SeqRecord):

    __slots__ = ('header', 'seq')

    def __init__(self, header : str, seq : str):
        self.header = header
        self.seq = seq
    # end def

    def __str__(self):
        seq_concise = self._get_consice_seq()
        return f'''header : {self.header},
seq: {seq_concise}.\n'''
    # end def

    def __repr__(self):
        seq_concise = self._get_consice_seq()
        return f'''Fasta(
    header={self.header!r},
    sequence={seq_concise!r}
)'''
    # end def

    def _get_consice_seq(self):
        n_chars_show = 30
        n_chars_omitted = len(self.seq) - 2*n_chars_show
        return '{}../{:,}bp/..{}'.format(
            self.seq[:n_chars_show],
            n_chars_omitted,
            self.seq[-n_chars_show:]
        )
    # end def
# end class
