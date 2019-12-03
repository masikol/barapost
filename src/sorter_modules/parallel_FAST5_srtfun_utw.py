# -*- coding: utf-8 -*-

def assign_version_2(fast5_list):
    # Assign version attribute to '2.0' -- multiFAST5
    for f5path in fast5_list:
        with h5py.File(f5path, 'a') as f5file:
            f5file.attrs["file_version"] = b"2.0"
        # end with
    # end for
# end def assign_version_2