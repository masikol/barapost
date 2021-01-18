# -*- coding: utf-8 -*-
# This module defines function that distributes input files among threads equally.

def spread_files_equally(fq_fa_list, n_thr):
    # Function distributes files among threads equally.
    # :param fq_fa_list: tuple (or list) of paths to files meant to be processed:
    # :type fq_fa_list: tuple<str>, list<str>;
    # :param n_thr: number of therads;
    # :type n_thr: int;

    # Create list of lists. It's length is equal to 'n_thr'
    multilist = [list() for _ in range(n_thr)]

    # Fill lists
    for i, fpath in enumerate(fq_fa_list):
        multilist[i % n_thr].append(fpath)
    # end for

    for sublist in multilist:
        yield sublist
    # end for
# end def spread_files_equally
