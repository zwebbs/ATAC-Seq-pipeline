# File Name: utils.py
# Created By: ZW
# Created On: 2022-09-30
# Purpose: Define utility functions for the ATAC-seq workflow.

# defined a function that returns 
def filter_seq_metadata(seq_list, target_key, **filter_kwargs):
    input_list = seq_list
    for key, value in filter_kwargs.items():
        out = list(filter(lambda i: i[key] == value, input_list))
        input_list = out
    return [o[target_key] for o in out]



