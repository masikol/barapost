# -*- coding: utf-8  -*-

def spread_files_equally(fq_fa_list, n_thr):
    """
    Function distributes files among processes equally.
    :param fq_fa_list: list of paths to files meant to be processed:
    :type fq_fa_list: list<str>;
    :param n_thr: number of therads to launch;
    :type n_thr: int;
    """

    glst = [list() for _ in range(n_thr)]

    for i, fpath in enumerate(fq_fa_list):
        glst[i % n_thr].append(fpath)
    # end for

    for sublist in glst:
        yield sublist
    # end for
# end def spread_files_equally