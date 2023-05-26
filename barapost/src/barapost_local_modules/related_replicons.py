# -*- coding: utf-8 -*-
# This module defines functions, which are necessary for searching for related replicons.

import re
import sys
from xml.etree import ElementTree

import socket
import http.client
from time import sleep

from src.platform import platf_depend_exit
from src.printlog import printlog_info, printlog_info_time, printlog_error, log_info
from src.printlog import printlog_warning, printlog_error_time, getwt
from src.lingering_https_get_request import lingering_https_get_request


def _get_record_title(record_id):
    # Function retrieves title (aka definition) and accession
    #   of a GenBank record by given accession or GI number.
    # :param record_id: accession or GI number of the record;
    # :type record_idi: str;
    # Returns tuple of two elements:
    #   (<RECORD_TITLE>, <RECORD_ACCESSION>)

    # We'll use E-utilities to communicate with GenBank
    eutils_server = "eutils.ncbi.nlm.nih.gov"
    esummary = "esummary.fcgi" # utility name

    # Configure URL
    url = "/entrez/eutils/{}?db=nuccore&id={}".format(esummary, record_id)

    # Sometimes (I never figured out why) this XML arrives empty, and StopIteration emerges.
    # So, if we just repeat this request, everything is going to be ok.
    error = True
    print_ok = False
    while error:
        # Send the request and get the response
        summary = lingering_https_get_request(eutils_server, url,
            "e-summary of nuccore record {}".format(record_id))

        # Parse XML that we've got
        root = ElementTree.fromstring(summary)

        # Elements of our insterest are all named "Item",
        #   but they have different tags.
        # They are children of element "DocSum", which is
        #   the first child of root
        try:
            docsum = next(iter(root))
        except StopIteration:
            print()
            printlog_info_time("Failed to retrieve data for record {}. Trying again...".format(record_id))
            print_ok = True # print this "ok" only after successful attepmt after fail
        else:
            if print_ok:
                printlog_info("ok")
            # end if
            error = False
        # end try
    # end while

    record_title = None
    record_acc = None

    # Search for title and accession
    for item in docsum.iter("Item"):
        if item.attrib["Name"] == "Title":
            record_title = item.text
        elif item.attrib["Name"] == "AccessionVersion":
            # Remove version just in case
            record_acc = re.search(r"(.*)\.[0-9]+", item.text).group(1)
        # end if
    # end for

    if record_title is None or record_acc is None:
        printlog_erro_time("Error 8989: can't access e-summary for `{}`".format(record_acc))
        platf_depend_exit(1)
    # end if

    return record_title, record_acc
# end get_record_title

class _NoIdentLabelError(Exception):
    pass
# end class NoIdentLabelError

class _NoLinkError(Exception):
    pass
# end class NoLinkError

class _NoAccError(Exception):
    pass
# end class NoAccError


def _is_redundant(nc_acc, accs):
    # Function checks if "NC-or-NW"-record is redundant (if it's non-RefSeq copy already exists in acc_dict).
    # :param nc_acc: accession number of NC-record;
    # :type nc_acc: str;
    # :param accs: tuple of accession numbers;
    # :type accs: tuple<str>;

    summary = lingering_https_get_request("www.ncbi.nlm.nih.gov",
        "/nuccore/{}?report=genbank&log$=seqview".format(nc_acc), "summary", nc_acc)

    try:
        # Find link to Identical GenBank Record

        # Firstly, get GI number of NC seqeunce:
        get_gi_url = "/nuccore/{}?report=gilist&log$=seqview&format=text".format(nc_acc)
        nc_gi_text = lingering_https_get_request("www.ncbi.nlm.nih.gov", get_gi_url,
            "GI of {}".format(nc_acc), nc_acc)
        nc_gi_text = nc_gi_text.replace('\n', '')
        nc_gi_re = re.search(r"\<pre\>([0-9]+).*\</pre\>", nc_gi_text)
        if nc_gi_re is None:
            raise _NoIdentLabelError("Error 771. Accession: {}. Please, contact the developer.".format(nc_acc))
        # end if

        nc_gi = nc_gi_re.group(1)

        # Retrieve identical GenBank sequence accession number.
        # NCBI redirects these requests and provides necessary location in headers.
        # So, we'll follow thin link.
        identical_gb_link = "/nuccore?LinkName=nuccore_nuccore_rsgb&from_uid={}".format(nc_gi)
        redirect_text = _ling_https_getreq_handl_301("www.ncbi.nlm.nih.gov", identical_gb_link,
            "link to identical genbank sequence", nc_acc)

        # Get accession number from the response text
        pattern = r"\<pre\>(.*).*\</pre\>"
        ident_acc_re = re.search(pattern, redirect_text.replace('\n', ''))

        if ident_acc_re is None:
            raise _NoIdentLabelError("Error 773. Accession: {}. Please, contact the developer.".format(nc_acc))
        # end if

        ident_acc = ident_acc_re.group(1).partition('.')[0]

    except (_NoIdentLabelError, _NoLinkError, _NoAccError) as err:
        printlog_error_time("Error: {}".format(err))
        platf_depend_exit(1)
    else:
        return ident_acc, ident_acc in accs
    # end try
# end def _is_redundant


class _DoesNotRedirectError(Exception):
    pass
# end class _DoesNotRedirectError


def _ling_https_getreq_handl_301(server, url, request_for=None, acc=None):
    # Name stands for "Lingering Https Get Request Handling 301".
    # Function performs a "lingering" HTTPS request.
    # It means that the function tries to get the response
    #     again and again if the request fails.
    # It handles 301-redirection in order to search for replicons related to "NC-records".
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
    while error:
        try:
            conn = http.client.HTTPSConnection(server, timeout=10) # create connection
            conn.request("GET", url) # ask for if there areresults
            response = conn.getresponse() # get the resonse

            # Handle redirection
            if response.code == 301:
                # Link to identical GenBank record is in "Location" header:
                redirect_url = response.getheader("Location")+"?report=accnlist&log$=seqview&format=text"
            else:
                raise _DoesNotRedirectError("NCBI does not redirect, although it must!")
            # end if

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
            printlog_warning("Can't connect to `{}`{}".format(server + url, comment_str))
            printlog_warning(str(err) )
            printlog_warning("the program will sleep for 30 seconds and try to connect again.")
            sleep(30)
        except _DoesNotRedirectError as err:
            printlog_error_time(str(err))
            printlog_error("Please, contact the developer.")
            platf_depend_exit(1)
        else:
            error = False # if no exception ocured, get out of the loop
        finally:
            conn.close()
        # end try
    # end while

    # And here goes simple "lingering_https_get_request",
    #   which will retrieve content from redirected location
    return lingering_https_get_request(server, redirect_url, request_for, acc)
# end def _ling_https_getreq_handl_301


def _deduplicate_replicons(repl_list, src_acc):
    # Function deduplicates related replicons: removes RefSeq-associated dublicates.
    # :param repl_list: list of discovered related replicons of format 'list<(ACC, DEFINITION)>';
    # :type repl_list: list<(str, str)>;
    # :param src_acc: 'source' accession, i.e. that from 'hits_to_download.tsv' of '-s' option;
    # :type src_acc: str;

    dedupl_list = list() # list of deduplicated pairs accession-definition
    accs = tuple(map(lambda x: x[0], repl_list))
    # redundant_lst = list() # collection of accession that will be removed from 'dedupl_list'

    # We will ntertain user -- show him/her this spinning thing (like conda does
    #   indicating that the script is actually working.
    krutiolka = ('|', '/', '-', '\\')
    krut_i = 0
    sys.stdout.write("\r {}".format(krutiolka[3]))
    sys.stdout.flush()

    for acc, hit_def in repl_list:

        if acc.startswith("NZ_") and acc[3:] in accs:
            # For "NZ"-records, check is simple: it's numerical part is identical to
            #   that from GenBank. Therefore we can just drop "NZ_" and go ahead.
            pass
        elif "NZ_" + acc == src_acc:
            # Simple non-NZ GenBank records should also be ignored here
            pass
        elif acc.startswith("NC_") or acc.startswith("NW_"):
            # "NC-or-NW"-records are no that simple:
            #   their numerical part is not identical to that from GenBank.
            # Therefore we need to check their equivalency by requesting NCBI site.
            # Example: this fungus: https://www.ncbi.nlm.nih.gov/nuccore/1800250711
            gb_acc, redundant = _is_redundant(acc, accs)
            if redundant and acc == src_acc:
                redund_tuples = filter(lambda x: x[0] == gb_acc, repl_list)
                rm_item = next(iter(redund_tuples))
                repl_list.remove(rm_item)
            # # end if
        elif acc != src_acc:
            # 'Source' accession are also in repl_list. Get rid of them.
            dedupl_list.append( (acc, hit_def) )
        # end if

        # Print this spinning thing
        sys.stdout.write("\r {}".format(krutiolka[krut_i]))
        sys.stdout.flush()
        krut_i = krut_i + 1 if krut_i != 3 else 0
    # end for

    # Print what we've got
    if len(dedupl_list) != 0:
        for i, (acc, hit_def) in enumerate(dedupl_list):
            print("\r  {}) {} - {}".format(i+1, acc, hit_def))
            log_info("  {}) {} - {}".format(i+1, acc, hit_def))
        # end for
    else:
        print("\r  It is the only replicon of this organism.")
        log_info("It is the only replicon of this organism.")
    # end if

    return dedupl_list
# end def _deduplicate_replicons


def _get_related_replicons(acc, acc_dict):
    # Generator finds replicons (other chromosomes or plasmids, sometimes even proviruses),
    #   which are related to Genbank record "discovered" by barapost-prober.py.
    #
    # :param acc: accession of a record "discovered" by barapost-prober.py;
    # :type acc: str;
    # :param acc_dict: dictionary {<ACCESSION>: <HIT_DEFINITION>};
    # :type acc_dict: dict<str: tuple<str>>;
    #
    # Yields tuples of a following structure:
    #     (<ACCESSION>, <RECORD_DEFINITION>)

    # We will save all titles in order not to duplicate records in our database
    repl_list = [(acc, acc_dict[acc])]

    # Elink utility returns links in DB_1, that are connected to given ID in DB_2
    eutils_server = "eutils.ncbi.nlm.nih.gov"
    elink = "elink.fcgi"

    # = Find BioSample ID =

    # Configure URL
    nuc2biosmp_url = "/entrez/eutils/{}?dbfrom=nuccore&db=biosample&id={}".format(elink, acc)

    # Get XML with our links
    text_link_to_bsmp = lingering_https_get_request(eutils_server,
        nuc2biosmp_url, "BioSample page", acc)

    # Parse this XML
    root = ElementTree.fromstring(text_link_to_bsmp)
    linkset = next(iter(root)).find("LinkSetDb")

    # XML should contain element "LinkSetDb"
    if linkset is None:
        printlog_warning("Cannot check replicons for `{}`: there is no BioSample page for this record.".format(acc))
        return list()
    # end if

    # Here we have BioSample ID
    biosmp_id = linkset.find("Link").find("Id").text

    # = Find assembly assotiated with this BioSample ID =

    # We will pass this BioSample ID through nuccore in order not to 
    #   allow requesting for over 7k transcripts, like for this fungus:
    #   https://www.ncbi.nlm.nih.gov/biosample/SAMN07457167
    # After this, only scaffolds (nearly 130 sequences) will be downloaded.

    # Configure URL
    biosmp2ass_url = "/entrez/eutils/{}?dbfrom=biosample&db=assembly&id={}".format(elink, biosmp_id)

    # Get XML with our links
    text_link_to_ass = lingering_https_get_request(eutils_server, biosmp2ass_url,
        "Assembly link assotiated with BioSample ID {}".format(biosmp_id))

    # Parse this XML
    root = ElementTree.fromstring(text_link_to_ass)
    linkset = next(iter(root)).find("LinkSetDb")

    # XML should contain element "LinkSetDb"
    if linkset is None:
        printlog_warning("""Cannot check replicons for `{}`: there is no assembly page for this record.""".format(acc))
        return list()
    # end if

    # Here we have BioSample ID
    ass_id = linkset.find("Link").find("Id").text

    # = Find GIs in nuccore assotiated with this Assembly ID =

    # Configure URL
    ass2nuc_url = "/entrez/eutils/{}?dbfrom=assembly&db=nuccore&id={}".format(elink, ass_id)

    # Get XML with our links
    text_link_to_nuc = lingering_https_get_request(eutils_server, ass2nuc_url,
        "Nucleotide links assotiated with assembly {}".format(ass_id))

    # Parse this XML
    root = ElementTree.fromstring(text_link_to_nuc)
    linkset = next(iter(root)).find("LinkSetDb")

    # XML should contain element "LinkSetDb"
    if linkset is None:
        printlog_error_time("""Cannot check replicons for `{}`:
  failed to find nuccore records for assembly {}.""".format(acc, ass_id))
        printlog_error("Please, contact the developer.")
        platf_depend_exit(1)
    # end if

    # We will ntertain user -- show him/her this spinning thing (like conda does
    #   indicating that the script is actually working.
    krutiolka = ('|', '/', '-', '\\')
    krut_i = 0
    sys.stdout.write("\r {}".format(krutiolka[3]))
    sys.stdout.flush()

    # Collect links
    for elem in linkset.iter():

        if elem.tag == "Id": # element "Id" contains our GI

            # Get GI, title and accession:
            rel_gi = elem.text
            rel_def, rel_acc = _get_record_title(rel_gi)

            # Print this spinning thing
            sys.stdout.write("\r {}".format(krutiolka[krut_i]))
            sys.stdout.flush()
            krut_i = krut_i + 1 if krut_i != 3 else 0

            # If accession is new -- update list
            if not rel_acc in map(lambda x: x[0], repl_list):
                # acc_dict[rel_acc] = rel_def # update acc_dict
                repl_list.append( (rel_acc, rel_def) )
            # end if
        # end if
    # end for
    return repl_list
# end def _get_related_replicons


def search_for_related_replicons(acc_dict):
    # Function searches for replicons related to those in 'hits_to_download.tsv'
    #   of specified with '-s' option.
    # :param acc_dict: dictionary comntaining accession data of hits;
    # :type acc_dict: dict<str: tuple<str, str, int>>;

    print()
    printlog_info_time("Searching for related replicons...")

    start_accs = tuple(acc_dict.keys()) # accessions, which were "discovered" by prober

    for i, acc in enumerate(start_accs):

        printlog_info("{}. {} ({}):".format(i+1, acc, acc_dict[acc]))

        # Search for related replicons:
        try:
            related_repls = _get_related_replicons(acc, acc_dict)
        except AttributeError:
            printlog_error_time("Parsing error: cannot find replicons related to {}.".format(acc))
            printlog_error("Please, contact the developer")
            platf_depend_exit(1)
        else:
            related_repls = _deduplicate_replicons(related_repls, acc)
        # end try
        for rel_acc, rel_def in related_repls:
            acc_dict[rel_acc] = rel_def
        # end for
    # end for

    print()
    if len(start_accs) != len(acc_dict): # there are some new replicons
        printlog_info_time("{} related replicons have been found.".\
            format(len(acc_dict) - len(start_accs)))
    else:
        printlog_info_time("No related replicons found.")
    # end if
    print()
# end def search_for_related_replicons
