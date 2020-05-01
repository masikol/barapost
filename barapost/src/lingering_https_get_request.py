# -*- coding: utf-8 -*-

import socket
import http.client
from time import sleep

from src.printlog import printl
from src.platform import platf_depend_exit


def lingering_https_get_request(server, url, logfile_path, request_for=None, acc=None):
    """
    Function performs a "lingering" HTTPS request.
    It means that the function tries to get the response
        again and again if the request fails.

    :param server: server address;
    :type server: str;
    :param url: the rest of url;
    :type url: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    :param request_for: some comment for error message;
    :type request_for: str;
    :param acc: GenBank accession;
    :type acc: str;

    Returns obtained response coded in UTF-8 ('str').
    """

    error = True
    while error:
        try:
            conn = http.client.HTTPSConnection(server, timeout=10) # create connection
            conn.request("GET", url) # ask for if there areresults
            response = conn.getresponse() # get the resonse

            if response.code != 200:
                printl(logfile_path, "Request failed with status code {}: {}".format(response.code, response.reason))
                platf_depend_exit(1)
            # end if

            resp_content = str(response.read(), "utf-8") # get response text
            conn.close()
        except (OSError, http.client.RemoteDisconnected, socket.gaierror, http.client.CannotSendRequest) as err:
            comment_str = ""
            if not request_for is None:
                comment_str += " requesting for {}".format(request_for)
                if not acc is None:
                    comment_str += " (accession: '{}')".format(acc)
                # end if
                comment_str += '.'
            # end if
            print()
            printl(logfile_path, "Can't connect to '{}'{}".format(server + url, comment_str))
            printl(logfile_path, str(err) )
            printl(logfile_path,"""the program will sleep for 30 seconds and try to connect again.""")
            sleep(30)
        else:
            error = False # if no exception ocured, get out of the loop
        # end try
    # end while
    return resp_content
# end def lingering_https_get_request