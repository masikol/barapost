
from src.Containers.SeqRecord import SeqRecord


class Fasta(SeqRecord):

    __slots__ = ['header', 'seq']

    def __init__(self, header : str, seq : str):
        self.header = header
        self.seq = seq
    # end def
# end class
