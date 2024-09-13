
from src.Containers.SeqRecord import SeqRecord


class Fasta(SeqRecord):

    __slots__ = ('header', 'seq')

    def __init__(self, header : str, seq : str):
        self.header = header
        self.seq = seq
    # end def

    def __str__(self):
        return f'''header : {self.header},
                seq: {self.seq}.\n'''
    # end def

    def __repr__(self):
        return f'Fasta(header={self.header!r}, sequence={self.seq!r})'
    # end def
# end class