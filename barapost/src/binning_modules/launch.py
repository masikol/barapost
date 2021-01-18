# -*- coding: utf-8 -*-
# Module defines functions that launch binning.

import multiprocessing as mp
from functools import partial

from src.spread_files_equally import spread_files_equally
from src.binning_modules.parallel_QA import init_paral_binning


def launch_single_thread_binning(fpath_list, str_func, tax_annot_res_dir, sens,
    min_qual, min_qlen, min_pident, min_coverage, no_trash):
    # Function launches single-thread binning, performed by finction 'srt_func'.
    #
    # :param fpath_list: list of path to files to process;
    # :type fpath_list: list<str>;
    # :param str_func: function that performs binning;
    # :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    # :type tax_annot_res_dir: str;
    # :param sens: binning sensitivity;
    # :type sens: str;
    # :param min_qual: threshold for quality filter;
    # :type min_qual: float;
    # :param min_qlen: threshold for length filter;
    # :type min_qlen: int (or None, if this filter is disabled);
    # :param min_pident: threshold for alignment identity filter;
    # :type min_pident: float (or None, if this filter is disabled);
    # :param min_coverage: threshold for alignment coverage filter;
    # :type min_coverage: float (or None, if this filter is disabled);
    # :param no_trash: loical value. True if user does NOT want to output trash files;
    # :type no_trash: bool;

    # Sort files in single thread:
    res_stats = map(
        partial(str_func, tax_annot_res_dir=tax_annot_res_dir, sens=sens,
            min_qual=min_qual, min_qlen=min_qlen, min_pident=min_pident,
            min_coverage=min_coverage, no_trash = no_trash),
        fpath_list)

    return res_stats
# end def launch_single_thread_binning


def launch_parallel_binning(fpath_list, str_func, tax_annot_res_dir, sens, n_thr,
    min_qual, min_qlen, min_pident, min_coverage, no_trash):
    # Function launches single-thread binning, performed by finction 'srt_func'.
    #
    # :param fpath_list: list of path to files to process;
    # :type fpath_list: list<str>;
    # :param str_func: function that performs binning;
    # :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    # :type tax_annot_res_dir: str;
    # :param sens: binning sensitivity;
    # :type sens: str;
    # :param n_thr: number of threads to launch;
    # :type n_thr: int;
    # :param min_qual: threshold for quality filter;
    # :type min_qual: float;
    # :param min_qlen: threshold for length filter;
    # :type min_qlen: int (or None, if this filter is disabled);
    # :param min_pident: threshold for alignment identity filter;
    # :type min_pident: float (or None, if this filter is disabled);
    # :param min_coverage: threshold for alignment coverage filter;
    # :type min_coverage: float (or None, if this filter is disabled);
    # :param no_trash: loical value. True if user does NOT want to output trash files;
    # :type no_trash: bool;

    # trick
    n_thr = min(n_thr, len(fpath_list))

    pool = mp.Pool(n_thr, initializer=init_paral_binning,
        initargs=(mp.Lock(), mp.Lock(),))

    res_stats = pool.starmap(partial(str_func,
            tax_annot_res_dir=tax_annot_res_dir,
            sens=sens,
            n_thr=n_thr,
            min_qual=min_qual,
            min_qlen=min_qlen,
            min_pident=min_pident,
            min_coverage=min_coverage,
            no_trash = no_trash),
        [(sublist,) for sublist in spread_files_equally(fpath_list, n_thr)])

    pool.close()
    pool.join()

    return res_stats
# end def launch_parallel_binning
