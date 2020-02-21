# -*- coding: utf-8 -*-
# This module defines functionsm which help to perform differently on different platforms.

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
    sys.exit(exit_code)
# end def platf_depend_exit

def get_logfile_path(script_name, outdir_path):
    """
    Function generates name of a log file
      according to name of "parent" program and time it was ran.

    :param script_name: name of program, which asks for log file path;
    :type script_name: str;
    :param outdir_path: path to output directory in which log file will be created;
    :type outdir_path: str;
    """

    # There are some troubles with file extention on Windows, so let's make a .txt file for it:
    log_ext = ".log" if not sys.platform.startswith("win") else ".txt"
    logfile_path = os.path.join(outdir_path, "{}_log_{}{}".format(script_name,
        strftime("%Y-%m-%d_%H-%M-%S", localtime(time())), log_ext))

    return logfile_path
# def get_logfile_path