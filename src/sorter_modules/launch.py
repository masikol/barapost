# -*- coding: utf-8 -*-

from functools import partial
import multiprocessing as mp

from src.sorter_modules.parallel_QA import init_paral_sorting
from src.spread_files_equally import spread_files_equally


def launch_single_thread_sorting(fpath_list, str_func, tax_annot_res_dir, sens,
    min_qual, min_qlen, logfile_path):

    # Sort files in single thread:
    res_stats = list(map(
        partial(str_func, tax_annot_res_dir=tax_annot_res_dir, sens=sens,
            min_qual=min_qual, min_qlen=min_qlen, logfile_path=logfile_path),
        fpath_list))

    return res_stats
# end def launch_single_thread_sorting


def launch_parallel_sorting(fpath_list, str_func, tax_annot_res_dir, sens, n_thr,
    min_qual, min_qlen, logfile_path):

    # trick
    n_thr = min(n_thr, len(fpath_list))

    print_lock = mp.Lock()
    write_lock = mp.Lock()

    pool = mp.Pool(n_thr, initializer=init_paral_sorting,
        initargs=(print_lock, write_lock))

    res_stats = pool.starmap(partial(str_func, tax_annot_res_dir=tax_annot_res_dir, sens=sens, n_thr=n_thr,
            min_qual=min_qual, min_qlen=min_qlen, logfile_path=logfile_path),
        [(sublist,) for sublist in spread_files_equally(fpath_list, n_thr)])

    pool.close()
    pool.join()

    return res_stats
# end def launch_parallel_sorting