# -*- coding: utf-8 -*-
# This module defines functions necessary for prober to interact with NCBI via net.

import os
from time import sleep
from re import search as re_search

import http.client
import urllib.parse

from src.lingering_https_get_request import lingering_https_get_request

from src.platform import platf_depend_exit
from src.printlog import printl, println, err_fmt, getwt, printn


def verify_taxids(taxid_list, logfile_path):
    """
    Funciton verifies TaxIDs passed to 'prober' with '-g' option.
    Function requests NCBI Taxonomy Browser and parses organism's name from HTML response.
    What is more, this function configures 'oraganisms' list - it will be included into BLAST submissions.

    :param taxid_list: list of TaxIDs. TaxIDs are strings, but they are verified to be integers
        during CL argument parsing;
    :type taxid_list: list<str>;
    :param logfile_path: path to logfile;
    :type logfile_path: str;

    Returns list of strings of the following format: "<tax_name> (taxid:<TaxID>)>"
    """
    organisms = list()
    if len(taxid_list) > 0:

        printl(logfile_path, "Verifying TaxIDs:")
        for taxid in taxid_list:
            println(logfile_path, "   {} - ".format(taxid))
            tax_url = "https://ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}".format(taxid)
            try:
                tax_resp = lingering_https_get_request("www.ncbi.nlm.nih.gov",
                    "/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}".format(taxid),
                    "taxonomy")
                tax_name = re_search(r"Taxonomy browser \((.+?)\)", tax_resp).group(1)
            except AttributeError:
                printl(logfile_path, "\aError: TaxID not found")
                print("Please, check your TaxID: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi")
                platf_depend_exit(1)
            except OSError as oserr:
                printl(logfile_path, err_fmt("something is wrong with connection:"))
                printl(logfile_path,  str(oserr) )
                platf_depend_exit(-2)
            else:
                printl(logfile_path, "\r   {} - {}".format(taxid, tax_name))
                organisms.append("{} (taxid:{})".format(tax_name, taxid))
            # end try
        # end for
        printl(logfile_path, '-'*30 + '\n')

    # end if
    return organisms
# end def verify taxids


def configure_request(packet, blast_algorithm, organisms, user_email):
    """
    Function configures the submissoin request to BLAST server.

    :param packet: FASTA_data_containing_query_sequences;
    :type packet: str;
    :param blast_algorithm: BLASTn algorithm to use;
    :type blast_algorithm: str;
    :param organisms: list of strings performing 'nt' database slices;
    :type organisms: list<str>;

    Returns a dict of the following structure:
    {
        "payload": the_payload_of_the_request (dict),
        "headers": headers of thee request (dict)
    }
    """

    payload = dict()
    payload["CMD"] = "PUT" # method
    payload["PROGRAM"] = "blastn" # program
    payload["MEGABLAST"] = "on" if "megablast" in blast_algorithm.lower() else "" # if megablast
    payload["BLAST_PROGRAMS"] = blast_algorithm # blastn algorithm
    payload["DATABASE"] = "nt" # db
    payload["QUERY"] = packet # FASTA data
    payload["HITLIST_SIZE"] = 1 # we need only the best hit
    if user_email != "":
        payload["email"] = user_email # user's email
        payload["tool"] = "barapost:_prober"
    # end if

    # 'nt' database slices:
    for i, org in enumerate(organisms):
        payload["EQ_MENU{}".format(i if i > 0 else "")] = org
    # end for

    payload["NUM_ORG"] = str( len(organisms) )
    payload = urllib.parse.urlencode(payload)
    headers = { "Content-Type" : "application/x-www-form-urlencoded" }

    return {"payload":payload, "headers": headers}
# end def configure_request


def send_request(request, pack_to_send, packet_size, packet_mode, filename, tmp_fpath, logfile_path):
    """
    Function sends a request to "blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        and then waits for satisfaction of the request and retrieves response text.

    :param request: request_data (it is a dict that 'configure_request()' function returns);
    :param request: dict<dict>;
    :param pack_to_send: current number (like id) of packet meant to be sent now.
    :type pack_to_send: int;
    :param pack_to_send: ordinal number of packet;
    :type pack_to_send: int;
    :param packet_size: numner of sequences in the packet;
    :type packet_size: int;
    :param logfile_path: path to logfile;
    :type logfile_path: str;

    Returns XML text of type 'str' with BLAST response.
    """
    payload = request["payload"]
    headers = request["headers"]

    server = "blast.ncbi.nlm.nih.gov"
    url = "/blast/Blast.cgi"
    error = True

    while error:
        try:
            conn = http.client.HTTPSConnection(server) # create a connection
            conn.request("POST", url, payload, headers) # send the request
            response = conn.getresponse() # get the response
            response_text = str(response.read(), "utf-8") # get response text
        except OSError as oserr:
            printl(logfile_path, "{} - 'https://blast.ncbi.nlm.nih.gov' is not available.".format(getwt()))
            printl(logfile_path,  str(oserr) )
            printl(logfile_path, "barapost will try to connect again in 30 seconds...\n")
            sleep(30)

        # if no exception occured
        else:
            error = False
        # end try
    # end while

    try:
        rid = re_search(r"RID = (.+)", response_text).group(1) # get Request ID
        rtoe = int(re_search(r"RTOE = ([0-9]+)", response_text).group(1)) # get time to wait provided by the NCBI server
    except AttributeError:
        printl(logfile_path, err_fmt("seems, ncbi has denied your request."))
        printl(logfile_path, "Response is in file 'request_denial_response.html'")
        with open("request_denial_response.html", 'w') as den_file:
            den_file.write(response_text)
        # end with
        platf_depend_exit(1)
    finally:
        conn.close()
    # end try

    # Save temporary data
    with open(tmp_fpath, 'w') as tmpfile:
        tmpfile.write("Request_ID: {}\n".format(rid))
        tmpfile.write("Packet_size: {}\n".format(packet_size))
        tmpfile.write("Packet mode: {}".format(packet_mode))
    # end with

    # Wait for results of alignment
    return( wait_for_align(rid, rtoe, pack_to_send, filename, logfile_path) )
# end def send_request


class BlastError:
    def __init__(self, code):
        self.code = code
    # end def
# end class BlastError


def wait_for_align(rid, rtoe, pack_to_send, filename, logfile_path):
    """
    Function waits untill BLAST server accomplishes the request.
    
    :param rid: Request ID to wait for;
    :type rid: str;
    :param rtoe: time in seconds estimated by BLAST server needed to accomplish the request;
    :type rtoe: int;
    :param pack_to_send: current packet (id) number to send;
    :type pack_to_send: int;
    :param filename: basename of current FASTA file;
    :type filename: str
    :param logfile_path: path to logfile;
    :type logfile_path: str;

    Returns XML response ('str').
    """

    printl(logfile_path, "\n{} - Requesting for current query status. Request ID: {},\n '{}'; Submission #{}".format(getwt(),
    rid, filename, pack_to_send[0],))
    # RTOE can be zero at the very beginning of resumption
    if rtoe > 0:

        printl(logfile_path, "{} - BLAST server estimates that alignment will be accomplished in {} seconds ".format(getwt(), rtoe))
        printl(logfile_path, "{} - Waiting for {}+3 (+3 extra) seconds...".format(getwt(), rtoe))
        # Server migth be wrong -- we will give it 3 extra seconds
        sleep(rtoe + 3)
        printl(logfile_path, "{} - {} seconds have passed. Checking if alignment is accomplished...".format(getwt(), rtoe+3))
    # end if

    server = "blast.ncbi.nlm.nih.gov"
    wait_url = "/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid

    whtspc_len = 6 + len("(requesting)")

    while True:
        resp_content = lingering_https_get_request(server, wait_url, logfile_path, "BLAST response")

        # if server asks to wait
        if "Status=WAITING" in resp_content:
            printn("\r{} - The request is being processed. Waiting{}{}".format(getwt(),
                ' '*whtspc_len, "\033[%dD" % whtspc_len))
            # indicate each 20 seconds with a dot
            for i in range(1, 7):
                sleep(10)
                printn("\r{} - The request is being processed. Waiting{}".format(getwt(), '.'*i))
            # end for
            printn("(requesting)")
            continue
        elif "Status=FAILED" in resp_content:
            # if job failed
            printl(logfile_path, '\n' + getwt() + " - Job failed\a\n")
            printl(logfile_path, "Resending this packet.")
            return None, BlastError(1)
        elif "Status=UNKNOWN" in resp_content:
            # if job expired
            printl(logfile_path, '\n' + getwt() + " - Job expired\a\n")
            printl(logfile_path, "Resending this packet.")
            return None, BlastError(1)
        # if results are ready
        elif "Status=READY" in resp_content:
            printl(logfile_path, "\n{} - Result for query '{}' #{} is ready!".format(getwt(), filename, pack_to_send[0]))
            # if there are hits
            if "ThereAreHits=yes" in resp_content:
                for i in range(15, 0, -5):
                    printl(logfile_path, '-' * i)
                # end for
                print("-\nRetrieving results...")

                # Retrieve human-readable text and put it into result directory
                retrieve_text_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&DESCRIPTIONS=1&ALIGNMENTS=1&RID=" + rid
                txt_align_res = lingering_https_get_request(server, retrieve_text_url, logfile_path,
                    "text version of BLAST response")

                # Count already existing plain text files in outdir:
                is_txt_response = lambda f: False if re_search(r"prober_blast_response_[0-9]+\.txt", f) is None else True
                outdir_path = os.path.dirname(logfile_path)
                response_num = len(tuple(filter(is_txt_response, os.listdir(outdir_path))))

                # Curent txt response file will have number 'response_num+1'
                txt_hpath = os.path.join(outdir_path, "prober_blast_response_{}.txt".format(response_num + 1))
                # Write text result for a human to read
                with open(txt_hpath, 'w') as txt_file:
                    txt_file.write(txt_align_res)
                # end with
            elif "ThereAreHits=no" in resp_content:
                # if there are no hits
                printl(logfile_path, getwt() + " - There are no hits. It happens.\n")
            else:
                # probably, job is failed if execution reaches here
                printl(logfile_path, '\n' + getwt() + " - Job failed\a\n")
                printl(logfile_path, "Resending this packet.")
                return None, BlastError(1)
            # end if
            break
        # end if
        # Execution should not reach here
        printl(logfile_path, '\n' + getwt() + " - Fatal error (-122). Please contact the developer.\a\n")
        platf_depend_exit(1)
    # end while

    # Retrieve XML result
    retrieve_xml_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&ALIGNMENTS=1&RID=" + rid
    xml_text = lingering_https_get_request(server, retrieve_xml_url, logfile_path,
        "XML BLAST response")

    if "Bad Gateway" in xml_text:
        printl(logfile_path, getwt() + " - Error! Bad Gateway! Data from last packet has been lost.")
        printl(logfile_path, "Resending this packet.")
        return None, BlastError(1)

    elif "Status=FAILED" in xml_text:
        printl(logfile_path, '\n' + getwt() + "BLAST error!: request failed")
        printl(logfile_path, "Resending this packet.")
        return None, BlastError(1)

    elif "to start it again" in xml_text:
        printl(logfile_path, '\n' + getwt() + "BLAST error!")
        printl(logfile_path, "Resending this packet.")
        return None, BlastError(1)

    elif "[blastsrv4.REAL]" in xml_text:
        printl(logfile_path, "BLAST server error:\n  {}".format(re_search(r"(\[blastsrv4\.REAL\].*\))", xml_text).group(1)))
        # Error code 2 indicated that we need to split packet and resubmit
        return None, BlastError(2)
    # end if

    return xml_text, BlastError(0)
# end def wait_for_align
