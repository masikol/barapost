from ..Containers.Fasta import Fasta

class FastaContainer():
    def __init__(self, file : Fasta, label : str):
        self.file = file
        self.label = label
    # end def
# end class

