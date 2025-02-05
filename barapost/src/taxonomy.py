# -*- coding: utf-8 -*-
# This module defines functions, which download and format organisms' taxonomys.

import re
import os
import glob

from src.lingering_https_get_request import lingering_https_get_request

from src.printlog import printlog_error, printlog_error_time
from src.platform import platf_depend_exit
from src.filesystem import remove_bad_chars
from src.filesystem import is_fasta, OPEN_FUNCS, FORMATTING_FUNCS, is_gzipped


ranks = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")

# These words at second (with index 1) position of title indicate that
#   actual species name are specified after it.
second_words_not_species = ("species", "sp.", "strain", "str.", "bacterium")

# Pattern for matching name of any rank except species:
high_tax_name_patt = r"[A-Z][a-z\.]+"
# Pattern to match species name
species_patt = r"[A-Za-z0-9\._]+"
# Pattern for matching whole taxonomy string.
# 6 semisolons with probable absence of name ending with species name.
# All without spaces.
proposed_fmt = r"(((%s)?;){6}(%s)?)" % (high_tax_name_patt, species_patt)

# Global list of accessions of sequences whose taxonomy is in taxonomy file
_tax_accs = list()


def init_tax_file(taxonomy_path):
    # Function for initializing taxonomy file (writing header to it)
    #
    # :param taxonomy_path: path to TSV taxonomy file;
    # :type taxonomy_path: str;

    with open(taxonomy_path, 'w') as tax_file:
        tax_file.write("{}\t{}\n".format("#ACCESSION",
            "TAXONOMY: {}".format(';'.join(ranks))))
    # end with
# end def init_tax_file


def fill_tax_accs(taxonomy_path):
    # Function reads TSV taxonomy file and adds all accessions stored in it
    #   to `_tax_accs`.
    #
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;

    global _tax_accs

    # Create this file if it does not exist
    if not os.path.exists(taxonomy_path):
        init_tax_file(taxonomy_path)
    # end if

    # Read the file, except the first line -- it's the header
    with open(taxonomy_path, 'r') as tax_file:
        tax_lines = tax_file.readlines()[1:]
    # end with

    # Extract the first column
    _tax_accs = list(map(lambda l: l.split('\t')[0], tax_lines))
# end def fill_tax_accs


def download_taxonomy(hit_acc, hit_def, taxonomy_path):
    # Function retrieves taxonomy of a hit from NCBI.
    # Moreover, it saves this taxonomy in file ``taxonomy_tsv:
    #     <accession>\t<taxonomy_str>
    #
    # :param hit_acc: hit accession;
    # :type hit_acc: str;
    # :param hit_def: definition of reference record;
    # :type hit_def: str;
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;

    # Get TaxID of the organism from GenBank summary:
    gb_summary = lingering_https_get_request("www.ncbi.nlm.nih.gov",
        "/nuccore/{}".format(hit_acc), "GenBank summary", hit_acc)


    re_search_obj = re.search(r"ORGANISM=([0-9]+)", gb_summary)
    if not re_search_obj is None:
        taxid = re_search_obj.group(1)

        # Get taxonomy page of the organism
        taxonomy_url = "/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}&lvl=3&lin=f&keep=1&srchmode=1&unlock".format(taxid)
        taxonomy_text = lingering_https_get_request("www.ncbi.nlm.nih.gov",
            taxonomy_url, "taxonomy", hit_acc)

        # This pattern will match taxonomic names along with their ranks
        tax_rank_pattern = r"TITLE=\"([a-z ]+)\"\>([A-Z].+?)\</a\>"

        # Get all taxonomic names of the organism
        taxonomy = re.findall(tax_rank_pattern, taxonomy_text)

        # We will convert ranks to lowercase just in case.
        # Firstly convert tuples to lists in order to change them:
        taxonomy = list(map(lambda x: list(x), taxonomy))

        # Remove odd information from beginnig of names:
        for i in range(len(taxonomy)):
            taxonomy[i][0] = taxonomy[i][0].lower() # just in case
        # end for

        # We will leave only following taxonomic ranks: domain, phylum, class, order, family, genus.
        # Species name requires special handling, it will be added later.
        ranks_to_select = ranks[:-1]

        # Remove redundant ranks:
        taxonomy = filter( lambda x: x[0].lower() in ranks_to_select, taxonomy )

        # Convert back to tuples:
        taxonomy = list(map(lambda x: tuple(x), taxonomy))

        # E.g., this record has no appropriate ranks: CP034535
        # Merely return it's definition
        if len(taxonomy) == 0:
            # Save taxonomy
            _tax_accs.append(hit_acc)
            with open(taxonomy_path, 'a') as tax_file:
                tax_file.write("{}\n".format('\t'.join( (hit_acc, hit_def) )))
            # end with
        # end if

        # Check if species name is specified like other ranks:
        check_direct_species_patt = r"TITLE=\"(species)\"\>([A-Za-z0-9 \.]+)\</a\>"
        match_direct_species = re.search(check_direct_species_patt, taxonomy_text)

        if not match_direct_species is None:
            # If species name is specified like other ranks, merely add it to list:
            taxonomy.append( (match_direct_species.group(1), match_direct_species.group(2).partition(" ")[2]) )
        else:
            # Otherwise we need to parse species name from title
            title = re.search(r"\<title\>Taxonomy browser \((.+)\)\</title\>", taxonomy_text).group(1)

            # Get words
            title = title.split(' ')

            # We will take all this words as species name.
            # Viruses also often have unpredictable names.
            #   Example: MN908947
            try:
                if title[1] in second_words_not_species or taxonomy[0][1].lower() == "viruses":
                    taxonomy.append( ("species", '_'.join(title[1:])) )
                else:
                    taxonomy.append( ("species", title[1]) )
                # end if
            except IndexError:
                # Handle absence of species name, e.g., this: AC150248.3
                # Well, nothing to append in this case!
                pass
            # end try
        # end if

        # Fill in missing ranks with empty strings
        for i in range(len(ranks)):
            if len(taxonomy) < i+1: # for this (missing in the end): AC150248
                taxonomy.append( (ranks[i], "") )
            elif taxonomy[i][0] != ranks[i]: # for this (mising in the middle): MN908947
                taxonomy.insert( i, (ranks[i], "") )
            # end if
        # end for

        # It will be a bit faster
        taxonomy = tuple(taxonomy)
    else:
        taxonomy = (
            ('Domain', '',),
            ('Phylum', '',),
            ('Class', '',),
            ('Order', '',),
            ('Family', '',),
            ('Genus', '',),
            ('Species', '',),
        )
    # end if

    # Save taxonomy
    _tax_accs.append(hit_acc)
    with open(taxonomy_path, 'a') as tax_file:
        tax_file.write("{}\n".format('\t'.join( (hit_acc, config_taxonomy_str(taxonomy)) )))
    # end with
# end def download_taxonomy


def find_taxonomy(hit_acc, hit_def, taxonomy_path):
    # Function returns taxonomy if it is already in taxonomy file
    #   and downloads it from NCBI Taxnomomy otherwise.
    #
    # :param hit_acc: hit accession;
    # :type hit_acc: str;
    # :param hit_def: definition of reference record;
    # :type hit_def: str;
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;

    if len(_tax_accs) == 0:
        fill_tax_accs(taxonomy_path)
    # end if

    # If hit is not new -- go further
    if hit_acc in _tax_accs:
        return
    # end if

    # If we've got a new accession -- download taxonomy
    download_taxonomy(hit_acc, hit_def, taxonomy_path)
# end def find_taxonomy


def parse_taxonomy(taxonomy_str):
    # Function parses ID of user's reference sequence and forms a taxonomy tuple
    #   if there is proper taxonomy string in ID line in fasta format.
    # "Proper" taxonomy string is following:
    #   '[ANYTHING BEFORE] <Domain>;<Phylum>;<Class>;<Order>;<Family>;<Genus>;<species> [ANYTHING AFTER]'
    # Spaces are not allowed. Ranks can be omitted in manner like this
    #   (order and species is missing):
    #   '[ANYTHING BEFORE] <Domain>;<Phylum>;<Class>;;<Family>;<Genus>; [ANYTHING AFTER]'
    # If there is no taxonomy string in sequence ID, we'll merely save this ID to taxonomy file.
    #
    # :param taxonomy_str: taxonomy string to parse;
    # :type taxonomy_str: str;
    #
    # Returns taxonomy tuple.

    # Check if `taxonomy_str` matches `proposed_fmt`
    proper_tax_match = re.search(proposed_fmt, taxonomy_str)

    # If there is a match and it taxonomic names are not empty,
    #   form taxonomic tuple:
    if not proper_tax_match is None and proper_tax_match.group(0) != ";"*(len(ranks)-1):

        tax_names = proper_tax_match.group(0).split(';')
        tax_names = tuple(map( str.strip, tax_names ))

        taxonomy = list()
        for i in range(len(ranks)):
            taxonomy.append( (ranks[i], tax_names[i]) )
        # end for
        taxonomy = tuple(taxonomy)
    # Otherwise we will merely use this sequence ID
    else:
        taxonomy = remove_bad_chars(taxonomy_str)
    # end if

    return taxonomy
# end def parse_taxonomy


def config_taxonomy_str(taxonomy_tuple):
    # Function converts taxonomy tuple (see function `parse_taxonomy`)
    #   and returns taxonomy string, where taxons are separated with `;`s
    # :param taxonomy_tuple: tuple which can be returned by function `parse_taxonomy`;
    # :type taxonomy_tuple: tuple<tuple<str, str>>;
    return ';'.join(map(lambda t: t[1], taxonomy_tuple))
# end def config_taxonomy_str


def get_tax_keys(taxonomy_path):
    # Function reads taxonomy file and returns it's
    #   content (just accessions, without taxonomy) as tuple.
    #
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;

    if not os.path.exists(taxonomy_path):
        init_tax_file(taxonomy_path)
    # end if

    # Read the file, except the first line - it's the header
    with open(taxonomy_path, 'r') as tax_file:
        tax_lines = tax_file.readlines()[1:]
    # end with

    # Extract the first column and return it as tuple
    tax_keys = tuple(map(lambda l: l.split('\t')[0], tax_lines))

    return tax_keys
# end def get_tax_dict


def get_tax_dict(taxonomy_path):
    # Function reads taxonomy file and returns it's
    #   content as dictionary.
    #
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;

    if not os.path.exists(taxonomy_path):
        init_tax_file(taxonomy_path)
    # end if

    tax_dict = dict()
    with open(taxonomy_path, 'r') as tax_file:
        tax_file.readline() # pass header
        for line in tax_file:
            acc, taxonomy = line.split('\t')
            tax_dict[acc] = parse_taxonomy(taxonomy.strip())
        # end for
    # end with

    return tax_dict
# end def get_tax_dict


def config_taxonomy_own_seq(taxonomy_str):
    # Function parses ID of user's reference sequence and forms a taxonomy tuple
    #   if there is proper taxonomy string in ID line in fasta format.
    # "Proper" taxonomy string is following:
    #   '[ANYTHING BEFORE] <Domain>;<Phylum>;<Class>;<Order>;<Family>;<Genus>;<species> [ANYTHING AFTER]'
    # Spaces are not allowed. Ranks can be omitted in manner like this
    #   (order and species is missing):
    #   '[ANYTHING BEFORE] <Domain>;<Phylum>;<Class>;;<Family>;<Genus>; [ANYTHING AFTER]'
    # If there is no taxonomy string in sequence ID, we'll merely save this ID to taxonomy file.
    #
    # :param taxonomy_str: taxonomy string to parse;
    # :type taxonomy_str: str;
    #
    # Returns taxonomy string.

    # Check if `taxonomy_str` matches `proposed_fmt`
    proper_tax_match = re.search(proposed_fmt, taxonomy_str)

    # If there is a match and it taxonomic names are not empty,
    #   form taxonomic tuple:
    if not proper_tax_match is None and proper_tax_match.group(0) != ";"*(len(ranks)-1):
        taxonomy = proper_tax_match.group(0)
    # Otherwise we will merely use this sequence ID
    else:
        taxonomy = remove_bad_chars(taxonomy_str)
    # end if

    return taxonomy
# end def config_taxonomy_own_seq


def save_taxonomy_directly(taxonomy_path, acc, taxonomy_str):
    # Function saves taxonomy to taxonomy file directly, without any downloading.
    #
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;
    # :param acc: accession of taxonomy entry;
    # :type acc: str;
    # :param taxonomy_str: taxonomy string to save;
    # :type taxonomy_str: is;

    if not os.path.exists(taxonomy_path):
        init_tax_file(taxonomy_path)
    # end if

    if len(_tax_accs) == 0:
        fill_tax_accs(taxonomy_path)
    # end if

    # Do not add redundant taxonomy line
    if not acc in _tax_accs:
        with open(taxonomy_path, 'a') as tax_file:
            tax_file.write("{}\n".format('\t'.join( (acc, config_taxonomy_own_seq(taxonomy_str)) )))
        # end with
        _tax_accs.append(acc)
    # end if
# end def def save_taxonomy_directly


def recover_taxonomy(acc, hit_def, taxonomy_path):
    # Function recovers missing taxonomy by given accession.
    #
    # :param acc: accession of taxonomy entry to recover;
    # :type acc: str;
    # :param hit_def: name of this sequence;
    # :type hit_def: sre;
    # :param taxonomy_path: path to TSV file with taxonomy;
    # :type taxonomy_path: str;

    if acc == "LAMBDA":
        # If we are missing lambda phage taxonomy -- just add it
        save_taxonomy_directly(taxonomy_path, acc, "Lambda-phage-nanopore-control")
    elif acc.startswith("OWN_SEQ_"):
        # If sequence is an "own seq" -- check fasta file

        # Get necessary title line from `local_seq_set.fasta`
        # Firstly find fasta file (it may be compressed)
        classif_dir = os.path.dirname(os.path.dirname(taxonomy_path))
        db_dir = os.path.join(classif_dir, "local_database")
        db_files = glob.glob("{}{}*".format(db_dir, os.sep))
        try:
            local_fasta = next(iter(filter(is_fasta, db_files)))
        except StopIteration:
            printlog_error_time("Error: cannot recover taxonomy for following sequence:")
            printlog_error(" `{} - {}`.".format(acc, hit_def))
            printlog_error("You can solve this problem by yourself (it's pretty simple).")
            printlog_error("Just add taxonomy line for {} to file `{}`".format(acc, taxonomy_path))
            printlog_error("  and run the program again.")
            platf_depend_exit(1)
        # end try

        # Find our line startingg with `acc`
        how_to_open = OPEN_FUNCS[is_gzipped(local_fasta)]
        fmt_func = FORMATTING_FUNCS[is_gzipped(local_fasta)]
        if is_gzipped(local_fasta):
            search_for = b">" + bytes(acc, 'ascii') + b" "
        else:
            search_for = ">{} ".format(acc)
        # end if

        with how_to_open(local_fasta) as fasta_file:
            for line in fasta_file:
                if line.startswith(search_for):
                    seq_name = fmt_func(line).partition(' ')[2] # get name of the sequence
                    save_taxonomy_directly(taxonomy_path, acc, seq_name)
                    break
                # end if
            # end for
        # end with
    else:
        # Try to find taxonomy in NCBI
        download_taxonomy(acc, hit_def, taxonomy_path)
    # end if
# end def check_taxonomy_consistensy
