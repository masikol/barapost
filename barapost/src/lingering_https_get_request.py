# -*- coding: utf-8 -*-

import socket
import http.client
from time import sleep

from src.printlog import printlog_info, printlog_error
from src.platform import platf_depend_exit

try:
    import ssl
except ImportError:
    pass
else:
    ssl._create_default_https_context = ssl._create_unverified_context
# end try


def lingering_https_get_request(server, url, request_for=None, acc=None):
    # Function performs a "lingering" HTTPS request.
    # It means that the function tries to get the response
    #     again and again if the request fails.
    #
    # :param server: server address;
    # :type server: str;
    # :param url: the rest of url;
    # :type url: str;
    # :param request_for: some comment for error message;
    # :type request_for: str;
    # :param acc: GenBank accession;
    # :type acc: str;
    #
    # Returns obtained response coded in UTF-8 ('str').

    error = True

    # We can get spurious 404 or sth due to instability of NCBI servers work.
    # Let's give it 3 attempts (with 15 sec spans in between),
    #   and if all them are unsuccessful -- teminate execution.
    attempt_i = 0
    max_attempts = 3

    while error:
        try:
            conn = http.client.HTTPSConnection(server, timeout=30) # create connection
            conn.request("GET", url) # ask for if there areresults
            response = conn.getresponse() # get the resonse

            if response.code != 200:
                if attempt_i < max_attempts and "ncbi.nlm.nih.gov" in server:
                    printlog_error("Error {}: {}.".format(response.code, response.reason))
                    printlog_error("It may be due to instable work of NCBI servers.")
                    printlog_error("{} attempts to connect left, waiting 15 sec..."\
                        .format(max_attempts - attempt_i))
                    attempt_i += 1
                else:
                    printlog_error("Cannot find {} for {}.".format(request_for, acc))
                    printlog_error("Request failed with status code {}: {}"\
                        .format(response.code, response.reason))
                    platf_depend_exit(1)
                # end if
            # end if

            resp_content = str(response.read(), "utf-8") # get response text
        except (OSError,\
                http.client.RemoteDisconnected,\
                socket.gaierror,\
                http.client.CannotSendRequest) as err:
            comment_str = ""
            if not request_for is None:
                comment_str += " requesting for {}".format(request_for)
                if not acc is None:
                    comment_str += " (accession: `{}`)".format(acc)
                # end if
                comment_str += '.'
            # end if
            print()
            printlog_info("Can't connect to `{}`{}".format(server + url, comment_str))
            printlog_info( str(err) )
            printlog_info("""the program will sleep for 30 seconds and try to connect again.""")
            sleep(30)
        else:
            error = False # if no exception ocured, get out of the loop
        finally:
            conn.close()
        # end try
    # end while
    return resp_content
# end def lingering_https_get_request
