# File Name: rparams.py
# Created By: ZW
# Created On: 2022-07-12
# Purpose: defines the dataclass holding rule parameters and resources


from dataclasses import dataclass
from typing import Dict, Union


class dotdict(dict):
        """dot.notation access to dictionary attributes"""
            __getattr__ = dict.get
                __setattr__ = dict.__setitem__
                    __delattr__ = dict.__delitem__




@dataclass(frozen=True)
class Rparam:
    rule_name: str
    parameters: Dict[str, Union[str,int, float, bool]]
    resources: Dict[str, Union[str, int]]



