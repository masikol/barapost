# -*- coding: utf-8 -*-
# This module defines functions, which download and format organisms' lineages.

import os
from time import sleep
from re import search as re_search, findall as re_findall

from src.lingering_https_get_request import lingering_https_get_request

from src.printlog import err_fmt
from src.platform import platf_depend_exit


ranks = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")

# These words at second (with index 1) position of title indicate that
#   actual species name are specified after it.
second_words_not_species = ("species", "sp.", "strain", "str.", "bacterium")

def remove_odd_info(some_name):
    """
    Function removes odd information from string meant to contain
      name of taxon higher than species.
    "Odd indormation" referrs to, e.g., word "Candidatus" for this P. ziziphi: CP025121.
    Fucntion makes an assumption that odd info is always writted at the beginning of string.

    :param some_name: string to be processed;
    :type some_name: str;
    """

    # If there are only one word in our name -- there are no odd info
    if some_name.count(' ') == 0:
        return some_name # return name itself
    # end if

    # Split into words
    names = some_name.split(' ')

    # Pattern matching taxon name higher than genus
    high_level_name_patt = r"^[A-Z][a-z]+$"
    # Function returns True if a word passed to it looks like high-level taxon name
    is_high_level_name = lambda w: not re_search(high_level_name_patt, w) is None

    # Get all high-level names in some_str
    high_level_names = tuple( filter(is_high_level_name, names) )

    if len(high_level_names) != 0:
        # Return the last one: all before it probably is odd
        return high_level_names[-1]
    else:
        # Well, we'll better merely return source string in this case
        return some_name
    # end if
# end def remove_odd_info


def download_lineage(hit_acc, hit_def, tax_file, logfile_path):
    """
    Function retrieves lineage of a hit from NCBI.
    Moreover, it saves this lineage in 'taxonomy' DBM file:
        {<accession>: <lineage_str>}

    :param hit_acc: hit accession;
    :type hit_acc: str;
    :param hit_def: definition of reference record;
    :type hit_def: str;
    :param tax_file: taxonomy file instance;
    :type tax_file: shelve.DbfilenameShelf;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    # Get TaxID of the organism from GenBank summary:
    gb_summary = lingering_https_get_request("www.ncbi.nlm.nih.gov",
        "/nuccore/{}".format(hit_acc), logfile_path, "GenBank summary", hit_acc)
    taxid = re_search(r"ORGANISM=([0-9]+)", gb_summary).group(1)

    # Get taxonomy page of the organism
    taxonomy_url = "/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}&lvl=3&lin=f&keep=1&srchmode=1&unlock".format(taxid)
    taxonomy_text = lingering_https_get_request("www.ncbi.nlm.nih.gov",
        taxonomy_url, logfile_path, "taxonomy", hit_acc)

    # This pattern will match taxonomic names along with their ranks
    tax_rank_pattern = r"TITLE=\"([a-z ]+)\"\>([A-Z].+?)\</a\>"

    # Get all taxonomic names of the organism
    lineage = re_findall(tax_rank_pattern, taxonomy_text)

    # We will convert ranks to lowercase just in case.
    # Firstly convert tuples to lists in order to change them:
    lineage = list(map(lambda x: list(x), lineage))


    # Remove odd information from beginnig of names:
    for i in range(len(lineage)):
        lineage[i][0] = lineage[i][0].lower() # just in case
        # Remove auxiliary odd data
        lineage[i][1] = remove_odd_info(lineage[i][1])
    # end for

    # We will leave only following taxonomic ranks.
    # Species name need special handling, it will be added later.
    ranks_to_select = ranks[:-1]

    # Remove redundant ranks:
    lineage = filter( lambda x: x[0].lower() in ranks_to_select, lineage )

    # Convert back to tuples:
    lineage = list(map(lambda x: tuple(x), lineage))

    # E.g., this record have no appropriate ranks: CP034535
    # Merely return it's definition
    if len(lineage) == 0:
        # Save taxonomy
        tax_file[hit_acc] = hit_def
        return hit_def
    # end if

    # Check if species name is specified like other ranks:
    check_direct_species_patt = r"TITLE=\"(species)\"\>([A-Za-z0-9 \.]+)\</a\>"
    match_direct_species = re_search(check_direct_species_patt, taxonomy_text)

    if not match_direct_species is None:
        # If species name is specified like other ranks, merely add it to list:
        lineage.append( (match_direct_species.group(1), match_direct_species.group(2).partition(" ")[2]) )
    else:
        # Otherwise we need to parse species name from title
        title = re_search(r"\<title\>Taxonomy browser \((.+)\)\</title\>", taxonomy_text).group(1)
        # Remove auxiliary odd data:
        genus_name = remove_odd_info(title)
        title = title.partition(genus_name)[1] + title.partition(genus_name)[2]

        # Get words
        title = title.split(' ')


        # We will take all this words as species name.
        # 
        # Viruses also often have unpredictable names.
        #   Example: MN908947
        try:
            if title[1] in second_words_not_species or lineage[0][1].lower() == "viruses":
                lineage.append( ("species", '_'.join(title[1:])) )
            else:
                lineage.append( ("species", title[1]) )
            # end if
        except IndexError:
            # Handle absence of species name, e.g., this: AC150248.3
            # Well, nothing to append in this case!
            pass
        # end try
    # end if

    # Fill in missing ranks with empty strings
    for i in range(len(ranks)):
        if len(lineage) < i+1: # for this (missing in the end): AC150248
            lineage.append( (ranks[i], "") )
        elif lineage[i][0] != ranks[i]: # for this (mising in the middle): MN908947
            lineage.insert( i, (ranks[i], "") )
        # end if
    # end for

    # It will be a bit faster
    lineage = tuple(lineage)

    # Save taxonomy
    tax_file[hit_acc] = lineage

    return get_str_to_print(lineage, hit_acc)
# end def download_lineage


def find_lineage(hit_acc, hit_def, tax_file, logfile_path):
    """
    Function returns lineage if it is already in taxonomy file
      and downloads it from NCBI Taxnomomy otherwise.

    :param hit_acc: hit accession;
    :type hit_acc: str;
    :param hit_def: definition of reference record;
    :type hit_def: str;
    :param tax_file: taxonomy file instance;
    :type tax_file: shelve.DbfilenameShelf;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    # Get all accessions in taxonomy file:
    tax_acc_exist = tuple(tax_file.keys())

    # If we've got a new accession -- download lineage
    if not hit_acc in tax_acc_exist:
        lineage = download_lineage(hit_acc, hit_def, tax_file, logfile_path)
    else:
        # If hit is not new -- simply retrieve it from taxonomy file
        lineage = get_str_to_print(tax_file[str(hit_acc)], hit_acc)
    # end if

    return lineage
# end def find_lineage


def get_str_to_print(lineage, hit_acc):
    """
    Funciton forms taxonomy string meant to be printed to conslole
      or written to classification file.
    Format of taxonomy string: "<Genus> <species>".
    Or if there are no genus and species ranks in taxonomy,
      it will be all present taxonomic ranks divided by semicolon.

    :param lineage: lineage object stored in taxonomy file;
    :type lineage: tuple<tuple<str>>;
    :param hit acc: accession of reference sequence;
    :type hit acc: str;
    """

    if isinstance(lineage, tuple):
        str_to_return = ""

        # Find genus and species names
        for rank, name in lineage:
            if rank == "genus":
                str_to_return = name + str_to_return # append to the beginniung
            elif rank == "species":
                str_to_return = str_to_return + " " + name # append to the end
            # end if
        # end for

        # If there are no genus and species, we'll return what we have
        if str_to_return == "" or str_to_return.startswith(' '):
            str_to_return = ';'.join(map(lambda x: x[1], lineage))
        # end if

        return str_to_return
    elif isinstance(lineage, str):
        return lineage
    else:
    # Execution must not reach here
        printl(logfile_path, "\nFatal error 8755.")
        printl(logfile_path, "Please, contact the developer -- it is his fault.")
        platf_depend_exit(1)
    # end if
# end def get_str_to_print


def get_lineage(hit_acc, tax_file, logfile_path):
    """
    Function returnes lineage by given accession from 'taxonomy' DBM file.

    :param hit_acc: hit accession;
    :type hit_acc: str;
    :param tax_file: taxonomy file instance;
    :type tax_file: shelve.DbfilenameShelf;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    # Retrieve taxonomy from taxonomy file
    try:
        if hit_acc in tax_file.keys():
            lineage = tax_file[hit_acc]
        # end if
    except KeyError:
        printl(logfile_path, err_fmt("{} is not in taxonomy file!".format(hit_acc)))
        printl(logfile_path, "Please, contact the developer.")
        platf_depend_exit(1)
    except OSError as oserr:
        printl(logfile_path, err_fmt(str(oserr)))
        platf_depend_exit(1)
    # end try

    # If we have beautiful formatted taxnonomy -- format it
    if isinstance(lineage, tuple):
        return get_str_to_print(lineage, hit_acc)
    # If we have unformatted custom sequence -- merely return this string
    elif isinstance(lineage, str):
        return lineage
    else:
    # Execution must not reach here
        printl(logfile_path, "\nFatal error 8753.")
        printl(logfile_path, "Please, contact the developer -- it is his fault.")
        platf_depend_exit(1)
    # end if
# end def get_lineage


def save_own_seq_taxonomy(seq_name, acc, tax_file):
    """
    Function parses ID of user's reference sequence and forms a taxonomy tuple
      if there is proper taxonomy string in ID line in fasta format.
    "Proper" taxonomy string is following:
      '[ANYTHING BEFORE] <Domain>;<Phylum>;<Class>;<Order>;<Family>;<Genus>;<species> [ANYTHING AFTER]'
    Spaces are not allowed. Ranks can be omitted in manner like this
      (order and species is missing):
      '[ANYTHING BEFORE] <Domain>;<Phylum>;<Class>;;<Family>;<Genus>; [ANYTHING AFTER]'
    If there is no taxonomy string in sequence ID, we'll merely save this ID to taxonomy file.

    :param seq_name: ID of reference sequencs;
    :type seq_name: str;
    :param acc: accession of reference sequence;
    :type acc: str;
    :param tax_file: taxonomy file instance;
    :type tax_file: shelve.DbfilenameShelf;
    """

    # Form complicated pattern to match taxonomy string

    # Pattern for matching name of any rank except species:
    high_tax_name_patt = r"[A-Z][a-z\.]+"
    # Pattern to match species name
    species_patt = r"[A-Za-z0-9\._]+"
    # Pattern for matching whole taxonomy string.
    # 6 semisolons with probable absence of name ending with species name.
    # All without spaces.
    proposed_fmt = r"(((%s)?;){6}(%s)?)" % (high_tax_name_patt, species_patt)

    # Find match
    proper_tax_match = re_search(proposed_fmt, seq_name)

    # If there is a match and it taxonomic names are not empty,
    #   form taxonomic tuple:
    if not proper_tax_match is None and proper_tax_match.group(0) != ";"*(len(ranks)-1):

        tax_names = proper_tax_match.group(0).split(';')
        tax_names = tuple(map( str.strip, tax_names ))

        lineage = list()
        for i in range(len(ranks)):
            lineage.append( (ranks[i], tax_names[i]) )
        # end for
        lineage = tuple(lineage)
    # Otherwise we will merely use this sequence ID
    else:
        lineage = seq_name.replace(' ', '_')
    # end if

    # Save taxonomy
    tax_file[acc] = lineage
# end def save_own_seq_taxonomy