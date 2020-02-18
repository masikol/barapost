# -*- coding: utf-8 -*-
# This module defines finctions, which are aimed to print information to console in different ways
#   duplicating this information in log file, if necessary.

import sys
from time import time, strftime, localtime, gmtime


start_time = time() # consider time of importing as start time

def getwt():
    """
    Function (get work time) returns time HH:MM:SS that has passed from start_time.
    """
    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def getwt


def get_full_time():
    """
    Function returns current time HH:MM:SS (YYYY-mm-dd HH:MM:SS).
    The first HH:MM:SS is from 'getwt'.
    """
    return getwt() + " ({}) ".format(strftime("%Y-%m-%d %H:%M:%S", localtime(time())))
# end def get_full_time


def printn(text=""):
    """
    Function prints text to the console without adding '\\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    """
    sys.stdout.write(text)
    sys.stdout.flush() # display text immediately
# end def printn


def printl(logfpath, text=""):
    """
    Function for printing text to console and to log file.
    """
    print(text)

    with open(logfpath, 'a') as logfile:
        logfile.write(str(text).strip('\r') + '\n')
        logfile.flush()
    # end with
# end def printl


def println(logfpath, text=""):
    """
    Function for printing text to console and to log file.
    The only difference from 'printl' -- text that is printed to console does not end with '\\n'
    """
    printn(text)

    with open(logfpath, 'a') as logfile:
        logfile.write(str(text).strip('\r') + '\n')
        logfile.flush()
    # end with
# end def printl


def err_fmt(text=""):
    """Function for configuring error messages"""
    return "\n  \a!! - ERROR: " + text + '\n'
# end def err_fmt
