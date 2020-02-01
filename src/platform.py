# -*- coding: utf-8 -*-

import sys
import os
from time import localtime, strftime, time

def platf_depend_exit(exit_code):
    """
    Function asks to press ENTER press on Windows
        and exits after that.

    :type exit_code: int;
    """
    if sys.platform.startswith("win"):
        input("Press ENTER to exit:")
    # end if
    exit(exit_code)
# end def platf_depend_exit

def get_logfile_path(script_name, outdir_path):

    # There some troubles with file extention on Windows, so let's make a .txt file for it:
    log_ext = ".log" if not sys.platform.startswith("win") else ".txt"
    logfile_path = os.path.join(outdir_path, "{}_log_{}{}".format(script_name,
        strftime("%Y-%m-%d_%H-%M-%S", localtime(time())), log_ext))

    return logfile_path
# def get_logfile_path