
from src.Containers.SeqRecord import SeqRecord

class Blow5(SeqRecord):
    def __init__(self, record : dict):
        self.record = dict()

        # Why do this?
        # Cause sametimes record['read_id'] has type == bytes (for some reason)
        # And pyslow5 raise exeption
        # This bullshit is an add hock, but problem sloved :)

        for key, value in record.items():
            if isinstance(key, bytes):
                new_key = key.decode('utf-8')
            else:
                new_key = key
            # end if
            self.record[new_key] = value
        # end for
    # end def
# end class

