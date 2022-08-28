"""" Show the parameters for a given model run.

This probably should be set up to be called from DataTablesVCN

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021, 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import pickle
import pprint
from collections import OrderedDict
from pathlib import Path
from typing import Union

import vcnmodel.model_params as model_params

MP = model_params.Params()
pp = pprint.PrettyPrinter()


def show_model_params(filename: Union[Path, str]):
    fn = Path(filename)
    with open(fn, "rb"):
        d = pickle.load(open(fn, "rb"))
    print(d.keys())
    if "Paramsx" in list(d.keys()):
        pp.pprint(d["Params"])

    else:
        # fill params from filename and modelpars
        print("\nmodelPars: ")
        pp.pprint(d["modelPars"])
        print()
        print("\nrunInfo: ")
        pp.pprint(d["runInfo"])
        print()
        print("\nDefault Params: ")
        pp.pprint(MP.Params)

    Par = MP.Params.copy()
    modelPars = d["modelPars"]
    runInfo = d["runInfo"]

    pk = list(Par.keys())
    for mpk in list(modelPars.keys()):
        if mpk in pk:
            if Par[mpk] == Par[mpk]:
                print(f"modelpar Duplicate key/value: ", mpk, Par[mpk], Par[mpk])
            else:
                print(
                    f"modelpar Matched key, different value: ", mpk, Par[mpk], Par[mpk]
                )
        else:
            print(f"Shoud add modelpar key {mpk:s} to Par keys")

    for mpk in list(runInfo.keys()):
        if mpk in pk:
            if Par[mpk] == Par[mpk]:
                print(f"runInfo Duplicate key/value: ", mpk, Par[mpk], Par[mpk])
            else:
                print(
                    f"runInfo Matched key, different value: ", mpk, Par[mpk], Par[mpk]
                )
        else:
            print(f"Shoud add runInfo key {mpk:s} to Par keys")
