# -*- coding: utf-8 -*-
# This module defines functions necessary for prober to interact with NCBI via net.

import os
import re
from time import sleep
import logging

import http.client
import urllib.parse

from src.lingering_https_get_request import lingering_https_get_request

from src.platform import platf_depend_exit
from src.printlog import log_info, printlog_info, printlog_info_time, printlog_error, printlog_error_time, getwt, printn


def verify_taxids(taxid_list):
    # Funciton verifies TaxIDs passed to prober with `-g` option.
    # Function requests NCBI Taxonomy Browser and parses organism's name from HTML response.
    # What is more, this function configures `oraganisms` list - it will be included into BLAST submissions.
    #
    # :param taxid_list: list of TaxIDs. TaxIDs are strings, but they are verified to be integers
    #     during CL argument parsing;
    # :type taxid_list: list<str>;
    #
    # Returns list of strings of the following format: "<tax_name> (taxid:<TaxID>)>"

    organisms = list()
    if len(taxid_list) > 0:

        printlog_info("Verifying TaxIDs:")
        for taxid in taxid_list:
            printn("   {} - ".format(taxid))
            try:
                tax_resp = lingering_https_get_request("www.ncbi.nlm.nih.gov",
                    "/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}".format(taxid),
                    "taxonomy")
                tax_name = re.search(r"Taxonomy browser \((.+?)\)", tax_resp).group(1)
            except AttributeError:
                printlog_error("\aError: TaxID not found")
                printlog_error("Please, check your TaxID: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi")
                platf_depend_exit(1)
            except OSError as oserr:
                printlog_error("Something is wrong with connection:")
                printlog_error( str(oserr) )
                platf_depend_exit(-2)
            else:
                print(tax_name)
                log_info("{} - {}".format(taxid, tax_name))
                organisms.append("{} (taxid:{})".format(tax_name, taxid))
            # end try
        # end for
        print('-'*30 + '\n')

    # end if
    return organisms
# end def verify taxids


def configure_request(packet, blast_algorithm, organisms, author_email):
    # Function configures the submissoin request to BLAST server.
    #
    # :param packet: FASTA_data_containing_query_sequences;
    # :type packet: str;
    # :param blast_algorithm: BLASTn algorithm to use;
    # :type blast_algorithm: str;
    # :param organisms: list of strings performing `nt` database slices;
    # :type organisms: list<str>;
    # :param author_email: author's email to send within request;
    # :type author_email: str;
    #
    # Returns a dict of the following structure:
    # {
    #     "payload": the_payload_of_the_request (dict),
    #     "headers": headers of thee request (dict)
    # }

    payload = dict()
    payload["CMD"] = "PUT" # method
    payload["PROGRAM"] = "blastn" # program
    payload["MEGABLAST"] = "on" if "megablast" in blast_algorithm.lower() else "" # if megablast
    payload["BLAST_PROGRAMS"] = blast_algorithm # blastn algorithm
    payload["DATABASE"] = "nt" # db
    payload["QUERY"] = packet # FASTA data
    payload["HITLIST_SIZE"] = 1 # we need only the best hit
    if author_email != "":
        payload["email"] = author_email # author's email
        payload["tool"] = "barapost_prober"
    # end if

    # `nt` database slices:
    for i, org in enumerate(organisms):
        payload["EQ_MENU{}".format(i if i > 0 else "")] = org
    # end for

    payload["NUM_ORG"] = str( len(organisms) )
    payload = urllib.parse.urlencode(payload)
    headers = { "Content-Type" : "application/x-www-form-urlencoded" }

    return {"payload":payload, "headers": headers}
# end def configure_request


def send_request(request, pack_to_send, packet_size, packet_mode, filename, tmp_fpath):
    # Function sends a request to "blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    #     and then waits for satisfaction of the request and retrieves response text.
    #
    # :param request: request_data (it is a dict that `configure_request()` function returns);
    # :param request: dict<dict>;
    # :param pack_to_send: current number (like id) of packet meant to be sent now.
    # :type pack_to_send: int;
    # :param pack_to_send: ordinal number of packet;
    # :type pack_to_send: int;
    # :param packet_size: numner of sequences in the packet;
    # :type packet_size: int;
    #
    # Returns XML text of type 'str' with BLAST response.

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
            printlog_info_time("`https://blast.ncbi.nlm.nih.gov` is not available.")
            printlog_info( str(oserr) )
            printlog_info("barapost will try to connect again in 30 seconds...\n")
            sleep(30)

        # if no exception occured
        else:
            error = False
        # end try
    # end while

    try:
        rid = re.search(r"RID = (.+)", response_text).group(1) # get Request ID
        rtoe = int(re.search(r"RTOE = ([0-9]+)", response_text).group(1)) # get time to wait provided by the NCBI server
    except AttributeError:
        printlog_error_time("Seems, NCBI has denied your request.")
        printlog_error("Response is in file `request_denial_response.html`")
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
        tmpfile.write("Packet_mode: {}".format(packet_mode))
    # end with

    # Wait for results of alignment
    return wait_for_align(rid, rtoe, pack_to_send, filename)
# end def send_request


class BlastError:
    def __init__(self, code):
        self.code = code
    # end def
# end class BlastError


def wait_for_align(rid, rtoe, pack_to_send, filename):
    # Function waits untill BLAST server accomplishes the request.
    #
    # :param rid: Request ID to wait for;
    # :type rid: str;
    # :param rtoe: time in seconds estimated by BLAST server needed to accomplish the request;
    # :type rtoe: int;
    # :param pack_to_send: current packet (id) number to send;
    # :type pack_to_send: int;
    # :param filename: basename of current FASTA file;
    # :type filename: str
    #
    # Returns XML response ('str').

    print()
    print("Requesting for current query status. Request ID: {}".format(rid))
    print(" `{}`; Submission #{}".format(filename, pack_to_send[0]))
    log_info("Requesting for current query status.")
    log_info("Request ID: {}; `{}`; Submission #{}".format(rid, filename, pack_to_send[0],))
    # RTOE can be zero at the very beginning of resumption
    if rtoe > 0:

        printlog_info_time("BLAST server estimates that alignment will be accomplished in {} seconds".format(rtoe))
        printlog_info_time("Waiting for {}+3 (+3 extra) seconds...".format(rtoe))
        # Server migth be wrong -- we will give it 3 extra seconds
        sleep(rtoe + 3)
        printlog_info_time("{} seconds have passed. Checking if alignment is accomplished...".format(rtoe+3))
    # end if

    server = "blast.ncbi.nlm.nih.gov"
    wait_url = "/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid

    whtspc_len = 6 + len("(requesting)")

    while True:
        resp_content = lingering_https_get_request(server, wait_url, "BLAST response")

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
            print()
            printlog_info_time("Job failed\a\n")
            printlog_info("Resending this packet.")
            return None, BlastError(2)
        elif "Status=UNKNOWN" in resp_content:
            # if job expired
            print()
            printlog_info_time("Job expired\a\n")
            printlog_info("Resending this packet.")
            return None, BlastError(1)
        # if results are ready
        elif "Status=READY" in resp_content:
            print()
            printlog_info("Result for query `{}` #{} is ready!".format(filename, pack_to_send[0]))
            # if there are hits
            if "ThereAreHits=yes" in resp_content:
                for i in range(15, 0, -5):
                    print('-' * i)
                # end for
                print("-\nRetrieving results...")

                # Retrieve human-readable text and put it into result directory
                retrieve_text_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&DESCRIPTIONS=1&ALIGNMENTS=1&RID=" + rid
                txt_align_res = lingering_https_get_request(server, retrieve_text_url, "text version of BLAST response")

                # Count already existing plain text files in outdir:
                is_txt_response = lambda f: not re.search(r"prober_blast_response_[0-9]+\.txt", f) is None
                outdir_path = os.path.dirname(logging.getLoggerClass().root.handlers[0].baseFilename) # tricky trick
                response_num = len(tuple(filter(is_txt_response, os.listdir(outdir_path))))

                # Curent txt response file will have number `response_num+1`
                txt_hpath = os.path.join(outdir_path, "prober_blast_response_{}.txt".format(response_num + 1))
                # Write text result for a human to read
                with open(txt_hpath, 'w') as txt_file:
                    txt_file.write(txt_align_res)
                # end with
            elif "ThereAreHits=no" in resp_content:
                # if there are no hits
                printlog_info_time("There are no hits. It happens.\n")
            else:
                # probably, job is failed if execution reaches here
                print()
                printlog_info_time("Job failed\a\n")
                printlog_info("Resending this packet.")
                return None, BlastError(2)
            # end if
            break
        # end if
        # Execution should not reach here
        printlog_error_time("Fatal error (-122). Please contact the developer.\a\n")
        platf_depend_exit(-122)
    # end while

    # Retrieve XML result
    retrieve_xml_url = "/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&ALIGNMENTS=1&RID=" + rid
    xml_text = lingering_https_get_request(server, retrieve_xml_url, "XML BLAST response")

    if "Bad Gateway" in xml_text:
        print()
        printlog_info_time("Bad Gateway. Data from last packet has been lost.")
        printlog_info("Resending this packet.")
        return None, BlastError(1)

    elif "Status=FAILED" in xml_text:
        print()
        printlog_info_time("BLAST error: request failed")
        printlog_info("Resending this packet.")
        return None, BlastError(2)

    elif "to start it again" in xml_text:
        print()
        printlog_info_time("BLAST error")
        printlog_info("Resending this packet.")
        return None, BlastError(2)

    elif "[blastsrv4.REAL]" in xml_text:
        blastsrv4_match = re.search(r"(\[blastsrv4\.REAL\].*\))", xml_text)
        blastsrv4_str = "" if blastsrv4_match is None else ": {}".format(blastsrv4_match.group(1))
        printlog_info_time("BLAST server error{}".format(blastsrv4_str))
        # Error code 2 indicated that we need to split packet and resubmit
        return None, BlastError(2)
    # end if

    return xml_text, BlastError(0)
# end def wait_for_align
