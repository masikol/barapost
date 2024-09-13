
from src.Containers.SeqRecord import SeqRecord

class Blow5(SeqRecord):

    __slots__ = ('record')

    def __init__(self, record : dict):
        self.record = record.copy()
    # end def
# end class

