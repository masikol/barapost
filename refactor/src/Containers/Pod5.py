
from src.Containers.SeqRecord import SeqRecord
from pod5 import ReadRecord

class Pod5(SeqRecord):

    __slots__ = ('record')

    def __init__(self, record : ReadRecord):
        self.record = record
    # end def
# end class

