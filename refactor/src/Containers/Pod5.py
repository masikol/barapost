
from src.Containers.SeqRecord import SeqRecord
from pod5 import ReadRecord


class Pod5(SeqRecord):

    __slots__ = ('record')

    def __init__(self, record : ReadRecord):
        self.record = record
    # end def

    def __str__(self):
        return f'record : {self.record}'
    # end def

    def __repr__(self):
        return f'''Pod5(
    record={self.record!r}
)'''
    # end def
# end class
