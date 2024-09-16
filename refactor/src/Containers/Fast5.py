
from h5py import File

from src.Containers.SeqRecord import SeqRecord

class Fast5(SeqRecord):

    __slots__ = ('out_file_handle', 'read_uuid')

    def __init__(self, 
                 out_file_handle : File, 
                 read_uuid : str):
        self.out_file_handle = out_file_handle
        self.read_uuid = read_uuid
    # end def

    def __str__(self):
        return f'''out_file_handle : {self.out_file_handle},
                read_uuid: {self.read_uuid}.\n'''
    # end def

    def __repr__(self):
        return f'Fast5(out_file_handle={self.out_file_handle!r}, read_uuid={self.read_uuid!r})'
    # end def
# end class