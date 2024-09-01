
from src.Containers.Fasta import Fasta
from src.Containers.ClassifContainer import ClassifContainer

# TODO: rename FastaContainer and FastQContainer
#   to sth like FastaClassificationContainter
class FastaContainer(ClassifContainer):
    def __init__(self, record : Fasta, label : str):
        self.record = record
        self.label = label
    # end def
# end class

