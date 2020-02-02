# -*- coding: utf-8 -*-

import os

def write_classification(res_tsv_lines, tsv_res_path):
    """
    Function writes result of blasting to result tsv file.

    :param res_tsv_lines: tsv lines returned by 'parse_align_results_xml()' funciton;
    :type res_tsv_lines: list<str>;
    :param tsv_res_path: path to reslut tsv file;
    :type tsv_res_path: str;
    """

    # If there is no result tsv file -- create it and write a head of the table.
    if not os.path.exists(tsv_res_path):
        with open(tsv_res_path, 'w') as tsv_res_file:
            tsv_res_file.write('\t'.join( ["QUERY_ID", "HIT_NAME", "HIT_ACCESSION", "QUERY_LENGTH",
                "ALIGNMENET_LENGTH", "IDENTITY", "GAPS", "E-VALUE", "AVG_PHRED33", "ACCURACY(%)"] ) + '\n')
        # end with
    # end if

    # Write reslut tsv lines to this file
    with open(tsv_res_path, 'a') as tsv_res_file:
        for line in res_tsv_lines:
            tsv_res_file.write(line + '\n')
        # end for
    # end with
# end def write_classification