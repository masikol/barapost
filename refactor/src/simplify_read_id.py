# This module defines function that simplifies ID of a sequence passed to it.
# This function handles Oxford Nanopore-like sequence IDs in a specific way.

import re

# According to
# https://github.com/nanoporetech/ont_h5_validator/blob/master/h5_validator/schemas/multi_read_fast5.yaml
_ont_read_pattern : re.Pattern = re.compile(
    r'([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})'
)


def simplify_read_id(read_id : str):
    srch_ont_read = re.search(_ont_read_pattern, read_id)
    if srch_ont_read is None:
        return read_id.partition(' ')[0][1:]
    else:
        return srch_ont_read.group(1)
# end def
