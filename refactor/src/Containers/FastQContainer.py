import math

from ..Containers.FastQ import FastQ

class FastQContainer():
    def __init__(self, file : FastQ, label : str, offset : int = 33):
        self.file = file
        self.label = label
        if offset not in [33, 64]:
            print(f"Unexpected offset = {offset}! Set offset in 33!")
            self.offset = 33
        else:
            self.offset = offset
    # end def

    def average_quality(self) -> float:
        total_pe = 0.0
        
        for char in self.file.quality:
            Q = ord(char) - self.offset
            pe = 10 ** (-Q / 10)
            total_pe += pe

        avg_pe = total_pe / len(self.file.quality)
        avg_Q = -10 * math.log10(avg_pe)
        return avg_Q
    # end def
# end class

