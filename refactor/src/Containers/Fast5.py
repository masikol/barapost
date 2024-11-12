
from h5py import File

from src.Containers.SeqRecord import SeqRecord


class Fast5(SeqRecord):

    __slots__ = ('file_handle', 'read_id')

    def __init__(self,
                 file_handle : File,
                 read_id : str):
        self.file_handle = file_handle
        self.read_id = read_id
    # end def

    def __str__(self):
        return f'''file_handle : {self.file_handle},
read_id: {self.read_id}.\n'''
    # end def

    def __repr__(self):
        return f'''Fast5(
    file_handle={self.file_handle!r},
    read_id={self.read_id!r}
)'''
    # end def
# end class
