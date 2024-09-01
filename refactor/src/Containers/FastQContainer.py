
from src.Containers.FastQ import FastQ
from src.Containers.ClassifContainer import ClassifContainer


class FastQContainer(ClassifContainer):
    def __init__(self, record : FastQ, label : str):
        self.record = record
        self.label = label
    # end def
# end class
