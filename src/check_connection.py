# -*- coding: utf-8 -*-

import urllib.request
from src.printlog import getwt, printn
from src.platform import platf_depend_exit

def check_connection(url):
    """
    Function checks if url is available.

    :return: None if url is available;
    """
    printn("Checking internet connection...")

    check_mark = "ok"

    try:
        status_code = urllib.request.urlopen(url).getcode()
        # Just in case
        if status_code != 200:
            print('\n' + getwt() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
            print("Check your internet connection.\a")
            print("Status code: {}".format(status_code))
            platf_depend_exit(-2)
        # end if
    except OSError as err:
        print('\n' + getwt() + " - Site 'https://blast.ncbi.nlm.nih.gov' is not available.")
        print("Check your internet connection.\a")
        print( str(err) )

        # 'urllib.request.HTTPError' can provide a user with information about the error
        if isinstance(err, urllib.request.HTTPError):
            print("Status code: {}".format(err.code))
            print(err.reason)
        # end if
        platf_depend_exit(-2)
    else:
        print("\rChecking internet connection... {}".format(check_mark))
    # end try
# end def check_connection