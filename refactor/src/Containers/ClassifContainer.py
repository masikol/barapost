
from src.Containers.SeqRecord import SeqRecord


class ClassifContainer:
    def __init__(self, record : SeqRecord, label : str):
        self.record = record
        self.label = label
    # end def
# end class
