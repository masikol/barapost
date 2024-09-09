
from h5py import File

from src.Containers.SeqRecord import SeqRecord

class Fast5(SeqRecord):

    __slots__ = ['out_file_handle', 'read_uuid']

    def __init__(self, 
                 out_file_handle : File, 
                 read_uuid : str):
        self.out_file_handle = out_file_handle
        self.read_uuid = read_uuid
    # end def
# end class

