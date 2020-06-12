# -*- coding: utf-8 -*-
# Module defines functions, which conifigure their return values according on filter options.

import os
from re import search as re_search
from src.filesystem import is_fasta


# |== Truth table ==|
#   for identity and coverage filters

# A = {Value is not minus};
# B = {Value passes filter};
# R = {Do NOT place sequence to trash file};
# 
# +---+---+---+
# | A | B | R |
# +===+===+===+
# | 0 | 0 | 1 |
# +---+---+---+
# | 0 | 1 | 1 |
# +---+---+---+
# | 1 | 0 | 0 |
# +---+---+---+
# | 1 | 1 | 1 |
# +---+---+---+
# It is just implication: R = (A => B) = (!A OR B)


def get_QL_filter(fpath, quality, length):
    """
    Function returns a filter function.
    This filter in turn returns True if annotation line passed to it passes all filters and False otherwise.
    Returns quality and length filter.

    :param fpath: path to input file;
    :type fpath: str;
    :param quality: threshold for quality filter;
    :type quality: float;
    :param length: threshold for length filter;
    :type length: int or None, if filter is disabled;
    """

    filters = list()

    # Ad quality filter
    # There will be minus instead of quality for fasta files
    if not is_fasta(fpath):

        # Someone can try to classify fasta files and then bin fast5 accoeding to the classification.
        # In this case TypeError will be raised. Then just return True.

        def qual_filter(x):
            try:
                return x[0] >= quality
            except TypeError:
                return True
            # end try
        # end def qual_filter

        filters.append(qual_filter)
    # end if

    # Add length filter
    if not length is None:
        filters.append(lambda x: x[1] >= length)
    # end if

    # Return "integral" filter
    return lambda line: all( (f(line) for f in filters) )
# end def get_QL_filter


def get_align_filter(pident, coverage):
    """
    Function returns a filter function.
    This filter in turn returns True if annotation line passed to it passes all filters and False otherwise.
    Returns alignment identity and coverage filter.

    :param pident: threshold for alignment identity filter;
    :type pident: float or None, if filter is disabled;
    :param coverage: threshold for alignment coverage filter;
    :type coverage: float or None, if filter is disabled;
    """

    filters = list()

    # Add identity filter
    # It will be placed to trash file only if identity is not "minus"
    #   and it passes filter. See Truth table above.
    if not pident is None:
        filters.append(lambda x: x[2] == '-' or x[2] / x[1] >= pident)
    # end if

    # Add coverage filter
    # It will be placed to trash file only if coverage is not "minus"
    #   and it passes filter. See Truth table above.
    if not coverage is None:
        filters.append(lambda x: x[3] == '-' or x[3] / x[1] >= coverage)
    # end if

    # Return "integral" filter
    return lambda line: all( (f(line) for f in filters) )
# end def get_align_filter

# Pattern will match .fasta, .fastq and .fast5 extentions without '.gz'
ext_pattern = r".*(\.(m)?f(ast)?(q|a|5))(\.gz)?$"


def get_QL_trash_fpath(fpath, outdir_path, quality, length):
    """
    Function configures path to trash (quality and length) file according to filter options.

    :param fpath: path to input file;
    :type fpath: str;
    :param outdir_path: path to output directory;
    :type outdir_path: str;
    :param quality: threshold for quality filter;
    :type quality: float;
    :param length: threshold for length filter;
    :type length: int or None, if filter is disabled;
    """

    trash_fpath = os.path.join(outdir_path, "trash")

    # Add quality to name name of trash file
    if not is_fasta(fpath):
        trash_fpath += "-q{}".format(quality)
    # end if

    # Add length to name name of trash file
    if not length is None:
        trash_fpath += "-m{}".format(length)
    # end if

    # Add appropriate extention
    trash_fpath += re_search(ext_pattern, fpath).group(1)

    return trash_fpath
# end def get_QL_trash_fpath


def get_align_trash_fpath(fpath, outdir_path, pident, coverage):
    """
    Function configures path to trash (identity and coverage) file according to filter options.

    :param fpath: path to input file;
    :type fpath: str;
    :param outdir_path: path to output directory;
    :param pident: threshold for alignment identity filter;
    :type pident: float or None, if filter is disabled;
    :param coverage: threshold for alignment coverage filter;
    :type coverage: float or None, if filter is disabled;
    """

    trash_fpath = os.path.join(outdir_path, "align_trash")

    # Add identity to name name of trash file
    if not pident is None:
        trash_fpath += "-i{}".format(round(pident*100,2))
    # end if

    # Add coverage to name name of trash file
    if not coverage is None:
        trash_fpath += "-c{}".format(round(coverage*100, 2))
    # end if

    # Add appropriate extention
    trash_fpath += re_search(ext_pattern, fpath).group(1)

    return trash_fpath
# end def get_align_trash_fpath
