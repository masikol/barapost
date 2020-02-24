# -*- coding: utf-8 -*-
# This module defines functions, which download and format organisms' lineages.

import os
import shelve
import urllib.request
from time import sleep
from re import search as re_search, findall as re_findall

from src.printlog import err_fmt
from src.platform import platf_depend_exit


ranks = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")


def download_lineage(hit_acc, taxonomy_path):
    """
    Function retrieves lineage of a hit from NCBI.
    Moreover, it saves this lineage in 'taxonomy' DBM file:
        {<accession>: <lineage_str>}

    :param hit_acc: hit accession;
    :type hit_acc: str;
    :param taxonomy_path: path to DBM file with taxonomy;
    :type taxonomy_path: str;
    """

    # Get TaxID of the organism from GenBank summary:
    gb_summary_url = "https://www.ncbi.nlm.nih.gov/nuccore/{}".format(hit_acc)
    gb_summary = urllib.request.urlopen(gb_summary_url).read().decode("utf-8")
    taxid = re_search(r"ORGANISM=([0-9]+)", gb_summary).group(1)

    # Get taxonomy page of the organism
    taxonomy_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}&lvl=3&lin=f&keep=1&srchmode=1&unlock".format(taxid)
    taxonomy_text = urllib.request.urlopen(taxonomy_url).read().decode("utf-8")

    # This pattern will match taxonomic names along with their ranks
    tax_rank_pattern = r"TITLE=\"([a-z ]+)\"\>([A-Z][a-z]+)\</a\>"

    # Get all taxonomic names of the organism
    lineage = re_findall(tax_rank_pattern, taxonomy_text)

    # We will convert ranks to lowercase just in case.
    # Firstly convert tuples to lists in order to change them:
    lineage = list(map(lambda x: list(x), lineage))
    for rank in lineage:
        rank[0] = rank[0].lower()
    # end for
    # Convert bck to tuples:
    lineage = list(map(lambda x: tuple(x), lineage))

    # We will leave only following taxonomic ranks.
    # Species name need special handling, it will be added later.
    ranks_to_select = ranks[:-1]

    # Remove redundant ranks:
    lineage = list(filter( lambda x: x[0].lower() in ranks_to_select, lineage ))

    # Check if species name is specified like other ranks:
    check_direct_species_patt = r"TITLE=\"(species)\"\>([A-Za-z0-9 \.]+)\</a\>"
    match_direct_species = re_search(check_direct_species_patt, taxonomy_text)

    if not match_direct_species is None:
        # If species name is specified like other ranks, merely add it to list:
        lineage.append( (match_direct_species.group(1), match_direct_species.group(2).partition(" ")[2]) )
    else:
        # Otherwise we need to parse species name from title
        title = re_search(r"\<title\>Taxonomy browser \((.+)\)\</title\>", taxonomy_text).group(1)
        title = title.split(' ')

        # These words at second (with index 1) indicate that
        #   actual species name are specified after it.
        second_words_not_species = ("species", "sp.", "strain", "str.")

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

    # It will be a bit faster
    lineage = tuple(lineage)

    # Save taxonomy
    with shelve.open(taxonomy_path, 'c') as tax_file:
        tax_file[hit_acc] = lineage
    # end with

    return get_str_to_print(lineage, hit_acc)
# end def get_lineage


def get_str_to_print(lineage, hit_acc):

    str_to_return = ""

    for rank, name in lineage:
        if rank == "genus":
            str_to_return = name + str_to_return
        elif rank == "species":
            str_to_return = str_to_return + " " + name
        # end if
    # end for

    # Presence of species name and absence of genus name is nonsence
    if str_to_return.startswith(' '):
        print("Taxonomy parsing error 1456")
        print("Please, contact the developer -- it is his fault.")
        print("Tell him the erroneous accession: '{}'".format(hit_acc))
        platf_depend_exit(1)
    # end if

    # If there are no genus and species, we'll return what we have
    if str_to_return == "":
        str_to_return = ';'.join(map(lambda x: x[1], lineage))
    # end if

    return str_to_return
# end def get_str_to_print


def find_lineage(hit_acc, taxonomy_path):
    """
    Function returns lineage if it is already in taxonomy file
      and downloads it from NCBI Taxnomomy otherwise.

    :param hit_acc: hit accession;
    :type hit_acc: str;
    :param taxonomy_path: path to DBM file with taxonomy;
    :type taxonomy_path: str;
    """

    # Get all accessions in taxonomy file:
    with shelve.open(taxonomy_path, 'c') as tax_file:
        tax_acc_exist = tuple(tax_file.keys())
    # end with

    # If we've got a new accession -- download lineage
    if not hit_acc in tax_acc_exist:
        lineage = download_lineage(hit_acc, taxonomy_path)
    else:
        # If hit is not new -- simply retrieve it from taxonomy file
        with shelve.open(taxonomy_path, 'r') as tax_file:
            lineage = get_str_to_print(tax_file[str(hit_acc)], hit_acc)
        # end with
    # end if

    return lineage
# end def find_lineage


def get_lineage(hit_acc, taxonomy_path):
    """
    Function returnes lineage by given accession from 'taxonomy' DBM file.

    :param hit_acc: hit accession;
    :type hit_acc: str;
    :param taxonomy_path: path to DBM file with taxonomy;
    :type taxonomy_path: str;
    """

    try:
        with shelve.open(taxonomy_path, 'r') as tax_file:
            if hit_acc in tax_file.keys():
                lineage = tax_file[hit_acc]
            # end if
        # end with
    except KeyError:
        print(err_fmt("{} is not in taxonomy file!".format(hit_acc)))
        print("Please, contact the developer.")
        platf_depend_exit(1)
    except OSError as oserr:
        print(err_fmt(str(oserr)))
        platf_depend_exit(1)
    # end try

    # If we have beautiful formatted taxnonomy
    if isinstance(lineage, tuple):
        return get_str_to_print(lineage, hit_acc)
    # If we have unformatted custom sequence
    elif isinstance(lineage, str):
        return lineage
    else:
    # Execution must not reach here
        print("\nFatal error 8753.")
        print("Please, contact the developer -- it is his fault.")
        platf_depend_exit(1)
    # end if
# end def get_lineage


def save_own_seq_taxonomy(seq_name, acc, taxonomy_path):

    high_tax_name_patt = r"[A-Z][a-z\.]+"
    species_patt = r"[A-Za-z0-9\._]+"
    proposed_fmt = r"((%s)?;){6}(%s)?" % (high_tax_name_patt, species_patt)

    proper_tax_match = re_search(proposed_fmt, seq_name)

    if not proper_tax_match is None:

        tax_names = proper_tax_match.group(1).split(';')
        tax_names = tuple(map( str.strip, tax_names ))

        lineage = list()
        for i in range(len(ranks)):
            if tax_names[i] != "":
                lineage.append( (ranks[i], tax_names[i]) )
            # end if
        # end for
    else:
        lineage = seq_name.replace(' ', '_')
    # end if

    with shelve.open(taxonomy_path, 'c') as tax_file:
        tax_file[acc] = lineage
    # end with
# end def save_own_seq_taxonomy