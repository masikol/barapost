# -*- coding: utf-8 -*-
# This module defines function that formats large number, dividing every 3 digits with underscore.

def get_undr_sep_number(number):
    """
    Function that formats large number, dividing every 3 digits with underscore.

    :param number: number to be formatted;
    :type numer: int;

    Returns formatted input number contverted to str.
    """
    undr_sep_num = str(number)
    for i in range(len(undr_sep_num)-4, -1, -4):
        undr_sep_num = undr_sep_num[: i+1] + '_' + undr_sep_num[i+1: ]
    # end for
    return undr_sep_num
# end def get_undr_sep_number