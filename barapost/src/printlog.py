# -*- coding: utf-8 -*-
# This module defines finctions, which are aimed to print information to console in different ways
#   duplicating this information in log file, if necessary.

import sys
import logging
from time import time, strftime, localtime, gmtime


start_time = time() # consider time of importing as start time

def getwt():
    # Function (get work time) returns time HH:MM:SS that has passed from start_time.

    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def getwt


def printlog_info_time(msg):
    print("{} - {}".format(getwt(), msg))
    logging.info("({} from start)\t{}".format(getwt(), msg))
# end def printlog

def printlog_info(msg):
    print(msg)
    logging.info("({} from start)\t{}".format(getwt(), msg))
# end def printlog


def log_info(msg):
    logging.info("({} from start)\t{}".format(getwt(), msg))
# end def printlog


def printlog_error_time(msg):
    print("\n{} - {}".format(getwt(), msg))
    logging.error("({} from start)\t{}.".format(getwt(), msg))
# end def printlog

def printlog_error(msg):
    print(msg)
    logging.error("({} from start)\t{}.".format(msg))
# end def printlog

def printlog_warning(msg):
    print(msg)
    logging.warning("({} from start)\t{}.".format(getwt(), msg))
# end def printlog


def get_full_time():
    # Function returns current time HH:MM:SS (YYYY-mm-dd HH:MM:SS).
    # The first HH:MM:SS is from 'getwt'.

    return getwt() + " ({})".format(strftime("%Y-%m-%d %H:%M:%S", localtime(time())))
# end def get_full_time


def printn(text=""):
    # Function prints text to the console without adding '\\n' in the end of the line.
    # Why not just to use 'print(text, end="")'?

    sys.stdout.write(text)
    sys.stdout.flush() # display text immediately
# end def printn
