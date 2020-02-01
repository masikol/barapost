# -*- coding: utf-8 -*-

import sys
from time import time, strftime, localtime, gmtime


start_time = time() # consider time of importing as start time

def getwt():
    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def getwt

def get_full_time():
    return getwt() + " ({}) ".format(strftime("%Y-%m-%d %H:%M:%S", localtime(time())))
# end def get_full_time

def printn(text=""):
    """
    Function prints text to the console without adding '\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
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
    return "\n   \a!! - ERROR: " + text + '\n'
# end def err_fmt
