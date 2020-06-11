# -*- coding: utf-8 -*-
# This module defines functions, which are necessary for searching for related replicons.

import re
from xml.etree import ElementTree

import socket
import http.client
from time import sleep

from src.platform import platf_depend_exit
from src.printlog import printl, printn, println, err_fmt, getwt
from src.lingering_https_get_request import lingering_https_get_request


def _get_record_title(record_id, logfile_path):
    """
    Function retrieves title (aka definition) and accession
      of a GenBank record by given accession or GI number.

    :param record_id: accession of GI number of the record;
    :type record_idi: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;

    Returns tuple of two elements:
      (<RECORD_TITLE>, <RECORD_ACCESSION>)
    """

    # We'll use E-utilities to communicate with GenBank

    eutils_server = "eutils.ncbi.nlm.nih.gov"
    esummary = "esummary.fcgi" # utility name

    # Configure URL
    url = "/entrez/eutils/{}?db=nuccore&id={}".format(esummary, record_id)

    # Send the request and get the response
    summary = lingering_https_get_request(eutils_server, url, logfile_path,
        "e-summary of nuccore record {}".format(record_id))

    # Parse XML that we've got
    root = ElementTree.fromstring(summary)

    # Elements of our insterest are all named "Item",
    #   but they have different tags.
    # They are children of element "DocSum", which is
    #   the first child of root
    docsum = next(iter(root.getchildren()))

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
        printl(logfile_path, "Error 8989: can't access e-summary for '{}'".format(record_acc))
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


def _is_redundant(nc_acc, accs, logfile_path):
    """
    Function checks if "NC"-record is redundant (if it's non'RefSeq copy already exists in acc_dict).

    :param nc_acc: accession number of NC-record;
    :type nc_acc: str;
    :param accs: tuple of accession numbers;
    :type accs: tuple<str>;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    summary = lingering_https_get_request("www.ncbi.nlm.nih.gov", "/nuccore/{}?report=genbank&log$=seqview".format(nc_acc),
        logfile_path, "summary", nc_acc)

    ident_label = "Identical GenBank Sequence"

    try:

        # Find link to Identical GenBank Record
        if not ident_label in summary:
            raise _NoIdentLabelError("Error 771. Accession: {}. Please, contact the developer.".format(nc_acc))
        # end if

        pattern = r"\<a class=\"brieflinkpopperctrl\" href=\"(/nuccore\?LinkName=nuccore_nuccore_rsgb&amp;from_uid=[0-9]+)\".*\>{}\</a\>".format(ident_label)
        link_re = re.search(pattern, summary)

        if link_re is None:
            raise _NoIdentLabelError("Error 772. Accession: {}. Please, contact the developer.".format(nc_acc))
        # end if
        link = link_re.group(1)

        # Retrieve identical GenBank sequence accession number.
        # NCBI redirects these requests and provides necessary location in headers.
        # So, we'll follow thin link.
        redirect_text = _ling_https_getreq_handl_301("www.ncbi.nlm.nih.gov", link,
            logfile_path, "link to {}".format(ident_label.lower()), nc_acc)

        # Get accession number from the response text
        pattern = r"\<pre\>(.*).*\</pre\>"
        ident_acc_re = re.search(pattern, redirect_text.replace('\n', ''))

        if ident_acc_re is None:
            raise _NoIdentLabelError("Error 773. Accession: {}. Please, contact the developer.".format(nc_acc))
        # end if

        ident_acc = ident_acc_re.group(1).partition('.')[0]

    except (_NoIdentLabelError, _NoLinkError, _NoAccError) as err:
        print('\n', str(err))
        platf_depend_exit(1)
    else:
        return ident_acc, ident_acc in accs
    # end try
# end def _is_redundant


class _DoesNotRedirectError(Exception):
    pass
# end class _DoesNotRedirectError


def _ling_https_getreq_handl_301(server, url, logfile_path, request_for=None, acc=None):
    """
    Name stands for "Lingering Https Get Request Handling 301".

    Function performs a "lingering" HTTPS request.
    It means that the function tries to get the response
        again and again if the request fails.
    It handles 301-redirection in order to search for replicons related to "NC-records".

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
            printl(logfile_path, "Can't connect to '{}'{}".format(server + url, comment_str))
            printl(logfile_path, str(err) )
            printl(logfile_path,"""the program will sleep for 30 seconds and try to connect again.""")
            sleep(30)
        except _DoesNotRedirectError as err:
            printl(logfile_path, str(err))
            printl(logfile_path, "Please, contact the developer.")
            platf_depend_exit(1)
        else:
            error = False # if no exception ocured, get out of the loop
        finally:
            conn.close()
        # end try
    # end while

    # And here goes simple "lingering_https_get_request",
    #   which will retrieve content from redirected location
    return lingering_https_get_request(server, redirect_url,
        logfile_path, request_for, acc)
# end def _ling_https_getreq_handl_301


def _deduplicate_replicons(repl_list, src_acc, logfile_path):
    """
    Function deduplicates related replicons: removes RefSeq-associated dublicates.
    :param repl_list: list of discovered related replicons of format 'list<(ACC, DEFINITION)>';
    :type repl_list: list<(str, str)>;
    :param src_acc: 'source' accession, i.e. that from 'hits_to_download.tsv' of '-s' option;
    :type src_acc: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    dedupl_list = list() # list of deduplicated pairs accession-definition
    accs = tuple(map(lambda x: x[0], repl_list))
    # redundant_lst = list() # collection of accession that will be removed from 'dedupl_list'

    for acc, hit_def in repl_list:

        if acc.startswith("NZ_") and acc[3:] in accs:
            # For "NZ"-records, check is simple: it's numerical part is identical to
            #   that from GenBank. Therefore we can just drop "NZ_" and go ahead.
            pass
        elif "NZ_" + acc == src_acc:
            # Simple non-NZ GenBank records should also be ignored here
            pass
        elif acc.startswith("NC_"):
            # "NC"-records are no that simple:
            #   their numerical part is not identical to that from GenBank.
            # Therefore we need to check their equivalency by requesting NCBI site.
            gb_acc, redundant = _is_redundant(acc, accs, logfile_path)
            if redundant and acc == src_acc:
                redund_tuples = filter(lambda x: x[0] == gb_acc, repl_list)
                rm_item = next(iter(redund_tuples))
                repl_list.remove(rm_item)
            # # end if
        elif acc != src_acc:
            # 'Source' accession are also in repl_list. Get rid of them.
            dedupl_list.append( (acc, hit_def) )
        # end if
    # end for

    # Print what we've got
    if len(dedupl_list) != 0:
        for i, (acc, hit_def) in enumerate(dedupl_list):
            printl(logfile_path, "  {}) {} - {}".format(i+1, acc, hit_def))
        # end for
    else:
        printl(logfile_path, "  It is the only replicon of this organism.")
    # end if

    return dedupl_list
# end def _deduplicate_replicons


def _get_related_replicons(acc, acc_dict, logfile_path):
    """
    Generator finds replicons (other chromosomes or plasmids, sometimes even proviruses),
      which are related to Genbank record "discovered" by barapost-prober.py.

    :param acc: accession of a record "discovered" by barapost-prober.py;
    :type acc: str;
    :param acc_dict: dictionary {<ACCESSION>: <HIT_DEFINITION>};
    :type acc_dict: dict<str: tuple<str>>;
    :param logfile_path: path to log file;
    :type logfile_path: str;

    Yields tuples of a following structure:
        (<ACCESSION>, <RECORD_DEFINITION>)
    """

    # We will save all titles in order not to duplicate records in our database
    repl_list = [(acc, acc_dict[acc])]

    # Elink utility returns links in DB_1, that are connected to given ID in DB_2
    eutils_server = "eutils.ncbi.nlm.nih.gov"
    elink = "elink.fcgi"

    # = Find BioSample ID =

    # Configure URL
    nuc2biosmp_url = "/entrez/eutils/{}?dbfrom=nuccore&db=biosample&id={}".format(elink, acc)

    # Get XML with our links
    text_link_to_bsmp = lingering_https_get_request(eutils_server, nuc2biosmp_url,
        logfile_path, "BioSample page", acc)

    # Parse this XML
    root = ElementTree.fromstring(text_link_to_bsmp)
    linkset = next(iter(root.getchildren())).find("LinkSetDb")

    # XML should contain element "LinkSetDb"
    if linkset is None:
        printl(logfile_path, "Cannot check replicons for '{}': \
there is no BioSample page for this record.".format(acc))
        return list()
    # end if

    # Here we have BioSample ID
    biosmp_id = linkset.find("Link").find("Id").text

    # = Find GIs in nuccore assotiated with this BioSample ID =

    # Configure URL
    biosmp2nuc_url = "/entrez/eutils/{}?dbfrom=biosample&db=nuccore&id={}".format(elink, biosmp_id)

    # Get XML with our links
    text_link_to_nuc = lingering_https_get_request(eutils_server, biosmp2nuc_url,
        logfile_path, "Nucleotide links assotiated with BioSample ID {}".format(biosmp_id))

    # Parse this XML
    root = ElementTree.fromstring(text_link_to_nuc)
    linkset = next(iter(root.getchildren())).find("LinkSetDb")

    # XML should contain element "LinkSetDb"
    if linkset is None:
        printl(logfile_path, """Cannot check replicons for '{}':
  there is no BioSample page for this record.""".format(acc))
        return list()
    # end if

    # Collect links
    for elem in linkset.iter():

        if elem.tag == "Id": # element "Id" contains our GI

            # Get GI, title and accession:
            rel_gi = elem.text
            rel_def, rel_acc = _get_record_title(rel_gi, logfile_path)

            # If accession is new -- update list
            if not rel_acc in map(lambda x: x[0], repl_list):
                # acc_dict[rel_acc] = rel_def # update acc_dict
                repl_list.append( (rel_acc, rel_def) )
            # end if
        # end if
    # end for
    return repl_list
# end def _get_related_replicons


def search_for_related_replicons(acc_dict, logfile_path):
    """
    Function searches for replicons related to those in 'hits_to_download.tsv'
      of specified with '-s' option.

    :param acc_dict: dictionary comntaining accession data of hits;
    :type acc_dict: dict<str: tuple<str, str, int>>;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    printl(logfile_path, "\n{} - Searching for related replicons...\n".format(getwt()))

    start_accs = tuple(acc_dict.keys()) # accessions, which were "discovered" by prober

    for i, acc in enumerate(start_accs):

        printl(logfile_path, "{}. {} ({}):".format(i+1, acc, acc_dict[acc]))

        # Search for related replicons:
        try:
            related_repls = _get_related_replicons(acc, acc_dict, logfile_path)
        except AttributeError:
            print("\nParsing error: cannot find replicons related to {}.".format(acc))
            print("Please, contact the developer")
            platf_depend_exit(1)
        else:
            related_repls = _deduplicate_replicons(related_repls, acc, logfile_path)
            for rel_acc, rel_def in related_repls:
                acc_dict[rel_acc] = rel_def
            # end for
        # end try
    # end for

    if len(start_accs) != len(acc_dict): # there are some new replicons
        printl(logfile_path, "\n{} - {} related replicons have been found.".format(getwt(),
            len(acc_dict) - len(start_accs)))
    else:
        printl(logfile_path, "\n{} - No related replicons found.\n".format(getwt()))
    # end if
# end def search_for_related_replicons