"""
This code is for "breaking the dendrites" from cells. Mostly it just allows us
to remove the dendrites from a hoc file. This is run as a script from the
command line. Run with python break_dendrites.py --help to see the command
arguments.

Note: *** DEPRECATED ***: replaced by compare_hoc.py (use that instead!)

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2020 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np
from neuron import h
from neuronvis import hoc_reader as HR
from neuronvis import swc_to_hoc

from vcnmodel.util.get_data_paths import get_data_paths

config = get_data_paths()

known_names = [
    "Axon_Hillock",
    "Axon_Initial_Segment",
    "Dendritic_Hub",
    "Dendritic_Swelling",
    "Distal_Dendrite",
    "Proximal_Dendrite",
    "soma",
]
basepath = config["cellDataDirectory"]

cells = [9]  # , 5, 6, 9, 10, 11, 13, 17, 30]

remove_types = ["Distal_Dendrite", "Dendritic_Swelling", "Dendritic_Hub"]
full_tag = "_Full"
nodend_tag = "_NoDend"
secmap = "sbem2"
# @dataclass
# class args:
#     pruneaxon: bool=False
#     prunedendrite: bool=True
#     prunedistal: bool=False
#     secmap: str='sbem'
#     topology: bool=True
#     filename: Union[str, Path] = ""

scales = {"x": 1.0, "y": 1.0, "z": 1.0, "r": 1.0, "soma": 1.0, "dend": 1.0}


def operate_on_cell(args):
    cells = [args.cell]
    for cell in cells:
        celln = f"VCN_c{cell:02d}"
        fn = Path(
            basepath, celln, "Morphology", Path(celln + full_tag).with_suffix(".swc")
        )
        if args.topology:
            ofile = None
        else:
            ofile = Path(
                basepath,
                celln,
                "Morphology",
                Path(celln + nodend_tag).with_suffix(".hoc"),
            )
        if fn.is_file():
            s = swc_to_hoc.SWC(filename=fn, secmap=secmap, scales=scales, args=args)
            hoc = s.write_hoc(filename=ofile)
            # print('hoc result: ', hoc)
            if args.topology:
                s.show_topology()

        # hoc = "\n".join(hoc)
        # r = h(hoc)
        # # h.topology()
        # # .show_topology()


def test_one():
    full_tag = "_Full"
    nodend_tag = "_NoDend"

    for cell in cells:
        celln = f"VCN_c{cell:02d}"
        fn = Path(
            basepath, celln, "Morphology", Path(celln + full_tag).with_suffix(".hoc")
        )

        ofile = Path(
            basepath, celln, "Morphology", Path(celln + nodend_tag).with_suffix(".hoc")
        )

        hoc = HR.HocReader(str(fn))
        # hoc.h.topology()
        sgr = hoc.retrieve_section_group()
        sl = hoc.get_section_lists()
        prefix = hoc.get_section_prefixes()

        secs_by_group = {key: [] for key in sl}

        for k, v in sgr.items():
            # print('v: ', v, k)
            secs_by_group[v].append(k)

        keep = {key: [] for key in sl if key not in remove_types}
        for k, v in secs_by_group.items():
            if k not in remove_types:
                for vx in v:
                    keep[k].append(vx)
            else:
                for fx in v:
                    sec = hoc.get_section(vx)
                # h.delete_section(sec)

        print(keep)


if __name__ == "__main__":
    print("This code is deprecated. Use compare_hoc.py instead")
    exit()
    parser = argparse.ArgumentParser(
        description="Manipulate and break dendrites or other parts off cell",
        argument_default=argparse.SUPPRESS,
        fromfile_prefix_chars="@",
    )
    parser.add_argument(
        dest="cell",
        type=int,
        default=None,
        help="Select the cell to convert (no default)",
    )
    parser.add_argument(
        "--prunedendrite",
        action="store_true",
        dest="prunedendrite",
        default=False,
        help="Prune all dendrite sections from the hoc output",
    )
    parser.add_argument(
        "--prunedistal",
        action="store_true",
        dest="prunedistal",
        default=False,
        help="Prune all dendrite sections proximal from the hoc output",
    )
    parser.add_argument(
        "--pruneaxon",
        action="store_true",
        dest="pruneaxon",
        default=False,
        help="Prune all axon sections from the hoc output",
    )
    parser.add_argument(
        "-t",
        "--topology",
        action="store_true",
        dest="topology",
        default=False,
        help="Show topology (blocks output)",
    )

    args = parser.parse_args()
    operate_on_cell(args)
