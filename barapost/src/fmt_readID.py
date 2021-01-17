# -*- coding: utf-8 -*-
# This module defines function that formats ID of a sequence passed to it.
# This function considers Oxford-Nanopore-like sequence IDs in a specific way.

import re

# According to
# https://github.com/nanoporetech/ont_h5_validator/blob/master/h5_validator/schemas/multi_read_fast5.yaml
ont_read_signature = r"([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})"

def fmt_read_id(read_id):

    srch_ont_read = re.search(ont_read_signature, read_id)
    if srch_ont_read is None:
        return '>' + read_id.partition(' ')[0][1:]
    else:
        return '>' + srch_ont_read.group(1)
# end def fmt_read_id