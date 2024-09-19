
from src.Containers.SeqRecord import SeqRecord

class Slow5(SeqRecord):

    __slots__ = ('record')

    def __init__(self, record : dict):
        self.record = record.copy()
    # end def

    def __str__(self):
        return f'record : {self.record}'
    # end def

    def __repr__(self):
        return f'''Slow5(
    record={self.record!r}
)'''
    # end def
# end class
