# File Name: utils.py
# Created By: ZW
# Created On: 2022-09-30
# Purpose: Define utility functions for the ATAC-seq workflow.

# defined a function that returns 
def filter_run_ids(seq_list, **kwargs):
    input_list = seq_list
    for key, value in kwargs.items():
        out = list(filter(lambda i: i[key] == value, input_list))
        input_list = out
    return [o["run_id"] for o in out]



