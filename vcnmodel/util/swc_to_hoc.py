"""
SWC File format from CNIC:

n T x y z R P n is an integer label that identifies the current point and
increments by one from one line to the next.

T is an integer representing the type of neuronal segment, such as soma, axon,
apical dendrite, etc. The standard accepted integer values are given below.

0 = undefined 1 = soma 2 = axon 3 = dendrite 4 = apical dendrite 5 = fork point
6 = end point 7 = custom x, y, z gives the cartesian coordinates of each node.

R is the radius at that node.

P indicates the parent (the integer label) of the current point or -1 to
indicate an origin (soma).

Python 3 version only 3-27-2019 pbm Handles Singleton "sections" in swc file by
inserting the last parent segment information.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import argparse
import datetime
from pathlib import Path

import numpy as np


# standard SWC types:
def swc_sectypes():
    swc_secs = {
        #  0: 'undefined',
        1: "soma",
        2: "axon",
        3: "basal_dendrite",
        4: "apical_dendrite",
        5: "custom",  # (user-defined preferences)
        6: "unspecified_neurites",
        7: "glia_processes",  # who knows why this is in here…
        10: "hillock",
        11: "unmyelinatedaxon",
        12: "dendrite",  # 'hub',
        # 13: 'proximal_dendrite',
        # 14: 'distal_dendrite',
    }
    return swc_secs


# section types for SBEM data on bushy cells (additional definitions)
def sbem_sectypes():
    sbem_sectype = {
        # new swc mapping
        0: "Undefined",
        1: "Soma",
        2: "Myelinated_Axon",
        3: "Basal_Dendrite",
        4: "Apical_Dendrite",
        5: "Custom",
        6: "Unspecified_Neurites",
        7: "Glia_Processe",
        8: "Blank",
        9: "Blank",
        10: "Axon_Hillock",
        11: "Unmyelinated_Axon",
        12: "Dendritic_Hub",
        13: "Proximal_Dendrite",
        14: "Distal_Dendrite",
        15: "Axon_Initial_Segment",
        16: "Axon_Heminode",
        17: "Axon_Node",
        18: "Dendritic_Swelling",
    }
    return sbem_sectype


# section types for SBEM data on bushy cells (additional definitions)
# This table is for swcs from Syglassfrom May 2021 (who changed it?)
def sbem2_sectypes():
    sbem2_sectype = {
        # new swc mapping
        0: "Undefined",
        1: "Soma",
        2: "Myelinated_Axon",
        3: "Basal_Dendrite",
        4: "Apical_Dendrite",
        5: "Custom",
        6: "Unspecified_Neurites",
        7: "Glia_Processe",
        8: "Blank",
        9: "Blank",
        10: "Axon_Hillock",
        11: "Dendritic_Swelling",
        12: "Dendritic_Hub",
        13: "Proximal_Dendrite",
        14: "Distal_Dendrite",
        15: "Axon_Initial_Segment",
        16: "Axon_Heminode",
        17: "Axon_Node",
    }
    return sbem2_sectype


# crenaming of cell parts to match cnmodel data tables (temporary)
renaming = {
    "basal_dendrite": "dendrite",
    "Basal_Dendrite": "dendrite",
    "Apical_Dendrite": "dendrite",
    "apical_dendrite": "dendrite",
    "proximal_dendrite": "dendrite",
    "Proximal_Dendrite": "dendrite",
    "distal_dendrite": "dendrite",
    "Distal_Dendrite": "dendrite",
    "Dendritic_Swelling": "dendrite",
    "hub": "dendrite",
    "Dendritic_Hub": "dendrite",
    "Axon_Hillock": "hillock",
    "Unmyelinated_Axon": "unmyelinatedaxon",
    "Axon_Initial_Segment": "initialsegment",
    "Axon_Heminode": "heminode",
    "Axon_Node": "node",
}


class SWC(object):
    """
    Encapsulates a morphology tree as defined by the SWC standard.

    Parameters
    ----------
    filename : str or None
        The name of an swc file to load
    types : dict or None
        A dictionary mapping {type_id: type_name} that describes the type IDs in
        the swc data (second column).
    data : ndarray or None
        Optionally, a data array may be provided instead of an swc file. This is
        used internally.
    scales : dict or None
        dict of format: {'x': 1.0, 'y': 1.0, 'z': 1.0, 'r': 1.0} to provide
        appropriate scaling along each of the axes.
    """

    def __init__(self, filename=None, types=None, secmap="swc", data=None, scales=None):
        self._dtype = [
            ("id", int),
            ("type", int),
            ("x", float),
            ("y", float),
            ("z", float),
            ("r", float),
            ("parent", int),
        ]

        self._id_lookup = None
        self._sections = None
        self._children = None
        self.scales = scales

        if secmap == "swc":
            self.sectypes = swc_sectypes
        elif secmap == "sbem":
            self.sectypes = sbem_sectypes()
        elif secmap == "sbem2":
            self.sectypes = sbem2_sectypes()
        else:
            raise ValueError("SWC number map type is not recognized: %s" % secmap)

        if types is not None:  # add-on or overwrite types to dictionary
            self.sectypes.update(types)

        if data is not None:
            self.data = data
        elif filename is not None:
            self.load(filename.with_suffix(".swc"))
            self.filename = filename
        else:
            raise TypeError("Must initialize with filename or data array.")

        self.sort()

    def load(self, filename):
        self.filename = filename
        print(f"Loading: {str(filename):s}")
        self.data = np.loadtxt(filename, dtype=self._dtype)
        if self.scales is not None:
            self.scale(
                x=self.scales["x"],
                y=self.scales["y"],
                z=self.scales["z"],
                r=self.scales["r"],
            )

    def copy(self):
        return SWC(data=self.data.copy(), types=self.sectypes)

    @property
    def lookup(self):
        """
        Return a dict that maps *id* to *index* in the data array.
        """
        if self._id_lookup is None:
            self._id_lookup = dict([(rec["id"], i) for i, rec in enumerate(self.data)])
            # self._id_lookup = {}
            # for i, rec in enumerate(self.data):
            # self._id_lookup[rec['id']] = i
        return self._id_lookup

    def children(self, ident):
        """
        Return a list of all children of the node *id*.
        """
        if self._children is None:  # build the child dict
            self._children = {}
            for rec in self.data:
                self._children.setdefault(rec["parent"], [])
                self._children[rec["parent"]].append(rec["id"])
        # print('children: ', self._children)
        return self._children.get(ident, [])

    def __getitem__(self, id):
        """
        Return record for node *id*.
        """
        return self.data[self.lookup[id]]

    def reparent(self, id):
        """
        Rearrange tree to make *id* the new root parent.
        """
        d = self.data

        # bail out if this is already the root
        if self[id]["parent"] == -1:
            return

        parent = -1
        while id != -1:
            oldparent = self[id]["parent"]
            self[id]["parent"] = parent
            parent = id
            id = oldparent

        self._children = None
        self.sort()

    @property
    def sections(self):
        """Return lists of IDs grouped by topological section.
        The first item in each list connects to the last item in a previous
        list.
        """
        # print('self.data: ', self.data)
        if self._sections is None:
            sections = []
            sec = []

            # find all nodes with nore than 1 child
            branchpts = set()
            endpoints = set(self.data["id"])
            endpoints.add(-1)
            seen = set()
            for r in self.data:
                p = r["parent"]
                if p in seen:
                    branchpts.add(p)
                else:
                    seen.add(p)
                    endpoints.remove(p)

            # build lists of unbranched node chains
            lasttype = self.data["type"][0]
            for r in self.data:
                sec.append(r["id"])
                if (
                    r["id"] in branchpts
                    or r["id"] in endpoints
                    or r["type"] != lasttype
                ):
                    sections.append(sec)
                    sec = []
                    lasttype = r["type"]

            self._sections = sections

        return self._sections

    def connect(self, parent_id, swc):
        """
        Combine this tree with another by attaching the root of *swc* as a
        child of *parent_id*.
        """
        data = swc.data.copy()
        shift = self.data["id"].max() + 1 - data["id"].min()
        data["id"] += shift
        rootmask = data["parent"] == -1
        data["parent"] += shift
        data["parent"][rootmask] = parent_id

        self.data = np.concatenate([self.data, data])
        self._children = None
        self.sort()

    def set_type(self, typ):
        self.data["type"] = typ

    def write_hoc(self, filename, types=None):
        """
        Write data to a HOC file.
        Each node type is written to a separate section list.
        """
        hoc = []
        # Add some header information
        hoc.extend([f"// Translated from SWC format by: swc_to_hoc.py"])
        hoc.append(f"// Source file: {str(self.filename):s}")
        hoc.append(f"// {datetime.datetime.now().strftime('%B %d %Y, %H:%M:%S'):s}")
        if self.scales is None:
            hoc.append(f"// No scaling")
        else:
            hoc.append(
                f"// Scaling: x: {self.scales['x']:f}, y: {self.scales['y']:f}, z: {self.scales['z']:f}, r: {self.scales['r']:f}"
            )
        hoc.append("")
        sectypes = self.sectypes.copy()
        print("sectypes: ", sectypes)
        for t in np.unique(self.data["type"]):
            print(t)
            if t not in sectypes:
                sectypes[t] = "type_%d" % t
        # create section lists
        screated = []
        for t in list(sectypes.values()):
            if t in screated:
                continue
            hoc.extend([f"objref {t:s}\n{t:s} = new SectionList()"])
            screated.append(t)
        hoc.append("")
        # create sections
        sects = self.sections

        hoc.append(f"create sections[{len(sects):d}]")
        sec_ids = {}

        for i, sec in enumerate(sects):
            # remember hoc index for this section
            endpt = self[sec[-1]]["id"]
            sec_id = len(sec_ids)
            sec_ids[endpt] = sec_id

            # add section to list
            hoc.append(f"access sections[{sec_id:d}]")
            typ = self[sec[0]]["type"]
            hoc.append(f"{sectypes[typ]:s}.append()")

            # connect section to parent
            p = self[sec[0]]["parent"]
            if p != -1:
                hoc.append(
                    f"connect sections[{sec_id:d}](0), sections[{sec_ids[p]:d}](1)"
                )

            # set up geometry for this section
            hoc.append("sections[%d] {" % sec_id)
            if len(sec) == 1:
                seg = sects[sec_ids[p]][-1]  # get last segement in the parent section
                rec = self[seg]
                if rec["r"] < 0.05:
                    print(f"MIN DIA ENCOUNTERED: {seg:d}, {rec['r']:f}")
                    rec["r"] = 0.05
                hoc.append(
                    f"  pt3dadd({rec['x']:f}, {rec['y']:f}, {rec['z']:f}, {rec['r']*2:f})  // seg={seg:d} Singleton repair: to section[sec_ids[p]:d]"
                )
            for seg in sects[sec_id]:
                rec = self[seg]
                if rec["r"] < 0.05:
                    print(f"MIN DIA ENCOUNTERED: {seg:d}, {rec['r']:f}")
                    rec["r"] = 0.05
                hoc.append(
                    f"  pt3dadd({rec['x']:f}, {rec['y']:f}, {rec['z']:f}, {rec['r']*2:f})   // seg={seg:d}"
                )
            hoc.append("}")

            hoc.append("")

        open(filename, "w").write("\n".join(hoc))
        print(f"Wrote {str(filename):s}")

    @property
    def root(self):
        """
        ID of the root node of the tree.
        """
        ind = np.argwhere(self.data["parent"] == -1)[0, 0]
        return self.data[ind]["id"]

    def sort(self):
        """
        Sort the tree in topological order.
        """
        order = self.branch(self.root)
        lt = self.lookup
        indexes = np.array([lt[i] for i in order], dtype=int)
        self.data = self.data[indexes]

        self._id_lookup = None
        self._sections = None

    def path(self, node):
        path = [node]
        while True:
            node = self[node]["parent"]
            if node < 0:
                return path
            path.append(node)

    def scale(self, x, y, z, r):
        self.data["x"] *= x
        self.data["y"] *= y
        self.data["z"] *= z
        self.data["r"] *= r

    def translate(self, x, y, z):
        self.data["x"] += x
        self.data["y"] += y
        self.data["z"] += z

    def branch(self, id):
        """
        Return a list of IDs in the branch beginning at *id*.
        """
        branch = [id]
        for ch in self.children(id):
            branch.extend(self.branch(ch))
        return branch

    def topology(self):
        """
        Print the tree topology.
        """
        path = []
        indent = ""
        secparents = [self[s[0]]["parent"] for s in self.sections]
        for i, sec in enumerate(self.sections):
            p = secparents[i]
            if p != -1:
                ind = path.index(p)
                path = path[: ind + 1]
                indent = indent[: (ind + 1) * 3]
            path.append(self[sec[-1]]["id"])

            # look ahead to see whether subsequent sections are children
            if p in secparents[i + 1 :]:
                this_indent = indent[:-2] + "├─ "
                indent = indent[:-2] + "│  │  "
            else:
                this_indent = indent[:-2] + "└─ "
                indent = indent[:-2] + "   │  "

            typ = self.sectypes[self[sec[0]]["type"]]
            if len(sec) > 10:
                secstr = "%s,...%s" % (
                    str(tuple(sec[:3]))[:-1],
                    str(tuple(sec[-3:]))[1:],
                )
            else:
                secstr = str(tuple(sec))
            print(
                "%ssections[%d] type=%s parent=%d %s" % (this_indent, i, typ, p, secstr)
            )


def main():
    # fn = 'VCN_Cells/VCN_c09/Morphology/VCN_c09_FullCell_edited.hoc'
    parser = argparse.ArgumentParser(
        description="swc_to_hoc file conversion",
        argument_default=argparse.SUPPRESS,
        fromfile_prefix_chars="@",
    )
    parser.add_argument(
        dest="input_file",
        action="store",
        default=None,
        help="Select theswc file to convert to hoc with section map",
    )
    parser.add_argument(
        "-s",
        "--secmap",
        type=str,
        default="swc",
        dest="secmap",
        choices=["swc", "sbem", "sbem2"],
        help="Choose section ampping)",
    )
    args = parser.parse_args()

    s = SWC(filename=Path(args.input_file), secmap=args.secmap)
    # s.topology()
    s.write_hoc(Path(args.input_file).with_suffix(".hocx"))


if __name__ == "__main__":
    main()
