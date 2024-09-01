
# TODO: remove

# from typing import MutableSequence

# from src.WriterSystem.FastaWriter import FastaWriter
# from src.WriterSystem.FastQWriter import FastQWriter
# from src.Containers.ClassifContainer import ClassifContainer


# class WriterSystem:

#     def __init__(self, _gzip_ : bool, n_max_out : int, _type_ : str):
#         self._gzip_ = _gzip_
#         self.n_max_out = n_max_out
#         self._type_ = _type_

#         if _type_ == 'fasta':
#             self.writer = FastaWriter(_gzip_, n_max_out)
#         elif _type_ == 'fastq':
#             self.writer = FastQWriter(_gzip_, n_max_out)
#         else:
#             raise ValueError(
#                 f'Invalid file type! Use fasta or fastq instead of {_type_}!'
#             )
#         # end if
#     # end def

#     def write(self, out_records : MutableSequence[ClassifContainer]):
#         self.writer.write(out_records)
#     # end def
# #end class
