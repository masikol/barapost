# -*- codingL utf-8 -*-

import os
import shelve
import urllib.request
from time import sleep
from re import search as re_search

from threading import Thread
import multiprocessing as mp

from src.filesystem import remove_tmp_files
from src.platform import platf_depend_exit
from src.printlog import err_fmt

own_seq_regex = r"OWN_SEQ_[0-9]+"


def download_lineage(gi, hit_def, acc, tax_annot_res_dir):
    """
    Function retrieves lineage of a hit from NCBI.
    It downloads INSDSeq XML file, since it is the smallest one among those containing lineage.
    Moreover, it saves this lineage in 'taxonomy' DBM file:
        {<accession>: <lineage_str>}

    :param gi: GI number of a hit;
    :type gi: str;
    :param hit_def: definition line of a hit;
    :type hit_def: str;
    :param acc: hit accession;
    :type acc: str;
    """

    indsxml_path = os.path.join(tax_annot_res_dir, "indsxml.gbc.xml")
    taxonomy_path = os.path.join(tax_annot_res_dir, "taxonomy","taxonomy")

    # Get all accessions in taxonomy file:
    tax_file = shelve.open(taxonomy_path, 'c')

    # If we've got a new accession -- download lineage
    if not re_search(own_seq_regex, acc) is None:
        tax_file[acc] = hit_def
    else:
        retrieve_url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbc_xml&id={}&".format(gi)

        def waiter(path):

            while not os.path.exists(path):
                sleep(0.1)
            # end while

            lin_regex = r"<INSDSeq_taxonomy>([A-Za-z0-9\.; ]+)</INSDSeq_taxonomy>"

            while re_search(lin_regex, open(path).read()) is None:
                sleep(0.1)
            # end while
        # end def waiter

        def downloader(retrieve_url, indsxml_path):
            error = True
            while error:
                try:
                    urllib.request.urlretrieve(retrieve_url, indsxml_path)
                except OSError:
                    print("\nError while requesting for lineage.\n Let's try again in 30 seconds.")
                    if os.path.exists(indsxml_path):
                        os.unlink(indsxml_path)
                    # end if
                    sleep(30)
                else:
                    error = False
                # end try
            # end while
        # end def downloader

        waiter_thread = Thread(target=waiter, args=(indsxml_path,))
        downloader_proc = mp.Process(target=downloader, args=(retrieve_url, indsxml_path))

        waiter_thread.start()
        downloader_proc.start()
        waiter_thread.join()
        downloader_proc.terminate()
        downloader_proc.join()

        # Get downloaded text:
        text = open(indsxml_path, 'r').read()

        try:

            # Find genus name and species name:
            org_name_regex = r"<INSDSeq_organism>([A-Z][a-z]+ [a-z]+(\. [a-zA-Z0-9]+)?)"
            org_name = re_search(org_name_regex, text).group(1)

            # Get species name:
            spec_name = re_search(r" ([a-z]+(\. [a-zA-Z0-9]+)?)", org_name).group(1)

            # Get full lineage
            lin_regex = r"<INSDSeq_taxonomy>([A-Za-z0-9;\. ]+)</INSDSeq_taxonomy>"
            lineage = re_search(lin_regex, text).group(1).strip('.')

            # Format of genus-species in lineage can vary.
            # We need to parse it correctly anyway.
            if spec_name != "sp": # no species info will be added to lineage for "Pseudarthrobacter sp."

                # Remove "Bacillus cereus group" from lineage
                for grp_cmplx in (" group", " complex"):
                    if grp_cmplx in lineage:
                        lineage = lineage[: lineage.rfind(';')]
                    # end if
                # end for

                # Check if genus and species names are in lineage (separated by space):
                gen_spec_regex = r"(([A-Z][a-z]+) [a-z]+(\. [a-zA-Z0-9]+)?).?$"
                gen_spec_in_lin = re_search(gen_spec_regex, lineage)

                # If there is genus and species in lineage
                if not gen_spec_in_lin is None:
                    genus_name = gen_spec_in_lin.group(2)
                    genus_count = lineage.count(" " + genus_name) # "...;Bacillus; Bacillus subtilis" are not allowed
                    if genus_count == 1:
                        pass # leave genus and species names intact
                    elif genus_count == 2: # remove odd genus name
                        lineage = lineage.replace(" " + genus_name + " ", " ")
                    else:
                        print(err_fmt("taxonomy parsing error"))
                        print("Please, contact the developer.")
                        print("Tell him that 'genus_count' is {}".format(genus_count))
                        platf_depend_exit(1)
                    # end if
                elif not ' ' + spec_name in lineage:
                    lineage += ";" + spec_name # if there are no species name -- add it
                # end if
            # end if

            # Remove all spaces:
            lineage = lineage.replace("sp. ", "sp._")
            lineage = lineage.replace("; ", ";")
            lineage = lineage.replace(" ", ";") # in orger to separate genus and species with semicolon

        except AttributeError:
            # If there is no correct lineage -- use hit definition instead of it

            # Format hit definition (get rid of stuff after comma)
            hit_def = hit_def[: hit_def.find(',')] if ',' in hit_def else hit_def
            hit_def = hit_def.replace(" complete genome", "") # sometimes there are no comma before it
            hit_def = hit_def.replace(' ', '_')

            lineage = hit_def
        # end try

        # Write lineage to taxonomy file
        tax_file[acc] = lineage
        remove_tmp_files(indsxml_path)
    # end if

    return lineage
# end def download_lineage


def get_lineage(hit_acc, taxonomy_path):
    """
    Function returnes lineage by given accession from 'taxonomy' DBM file.

    :param hit_acc: hit accession;
    :type hit_acc: str;
    """

    try:
        with shelve.open(taxonomy_path, 'r') as tax_file:
            if hit_acc in tax_file.keys():
                lineage = tax_file[hit_acc]
                return lineage
            else:
                # 'lcl|ACCESSION...' entries can be given with '.1'
                #   (or '.2', whatever) terminus by blastn.
                # There is no '.1' terminus in taxonomy file.
                # If we prune this terminus, key can hit the annotation.
                hit_acc_pruned = re_search(r"^(.+)\.[0-9]+$", hit_acc)
                if not hit_acc_pruned is None:
                    hit_acc_pruned = hit_acc_pruned.group(1)
                    if hit_acc_pruned in tax_file.keys():
                        lineage = tax_file[hit_acc_pruned]
                        return lineage
                    # end if
                # end if
                raise KeyError
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
# end def get_lineage