# -*- coding: utf-8 -*-
# Module defines functions that launch sorting.

from functools import partial
import multiprocessing as mp

from src.sorter_modules.parallel_QA import init_paral_sorting
from src.spread_files_equally import spread_files_equally


def launch_single_thread_sorting(fpath_list, str_func, tax_annot_res_dir, sens,
    min_qual, min_qlen, logfile_path):
    """
    Function launches single-thread sorting, performed by finction 'srt_func'.

    :param fpath_list: list of path to files to process;
    :type fpath_list: list<str>;
    :param str_func: function that performs sorting;
    :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    :type tax_annot_res_dir: str;
    :param sens: sorting sensitivity;
    :type sens: str;
    :param min_qual: minimum quality to keep;
    :type min_qual: float;
    :param min_qlen: minimmum sequence length to keep;
    :type min_qlen: int (or None, if this feature is disabled);
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    # Sort files in single thread:
    res_stats = map(
        partial(str_func, tax_annot_res_dir=tax_annot_res_dir, sens=sens,
            min_qual=min_qual, min_qlen=min_qlen, logfile_path=logfile_path),
        fpath_list)

    return res_stats
# end def launch_single_thread_sorting


def launch_parallel_sorting(fpath_list, str_func, tax_annot_res_dir, sens, n_thr,
    min_qual, min_qlen, logfile_path):
    """
    Function launches single-thread sorting, performed by finction 'srt_func'.

    :param fpath_list: list of path to files to process;
    :type fpath_list: list<str>;
    :param str_func: function that performs sorting;
    :param tax_annot_res_dir: path to directory containing taxonomic annotation;
    :type tax_annot_res_dir: str;
    :param sens: sorting sensitivity;
    :type sens: str;
    :param n_thr: number of threads to launch;
    :type n_thr: int;
    :param min_qual: minimum quality to keep;
    :type min_qual: float;
    :param min_qlen: minimmum sequence length to keep;
    :type min_qlen: int (or None, if this feature is disabled);
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    # trick
    n_thr = min(n_thr, len(fpath_list))

    print_lock = mp.Lock() # lock for printing to console
    write_lock = mp.Lock() # lock for writing to result file(s)

    pool = mp.Pool(n_thr, initializer=init_paral_sorting,
        initargs=(print_lock, write_lock))

    res_stats = pool.starmap(partial(str_func,
        tax_annot_res_dir=tax_annot_res_dir,
        sens=sens,
        n_thr=n_thr,
        min_qual=min_qual,
        min_qlen=min_qlen,
        logfile_path=logfile_path), [(sublist,) for sublist in spread_files_equally(fpath_list, n_thr)])

    pool.close()
    pool.join()

    return res_stats
# end def launch_parallel_sorting