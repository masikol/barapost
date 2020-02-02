# -*- coding: utf-8 -*-

def get_undr_sep_number(number):
    undr_sep_num = str(number)
    for i in range(len(undr_sep_num)-4, -1, -4):
        undr_sep_num = undr_sep_num[: i+1] + '_' + undr_sep_num[i+1: ]
    # end for
    return undr_sep_num
# end def get_undr_sep_number