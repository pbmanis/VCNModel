"""
This module takes one hoc file description of a neuron, and compares it to
another one. The file written out is trimmed version of the first file, with the
dendrite sections that are missing in the second file removed. Elements such as
the soma and various axon parts are kept intact. In addition, the keepids list
(see below) identifies IDs that need to be kept for the output file to be intact
(they might be missing in the comparison file).

This module is used to remove certain dendrites (un-innervated) from a hoc file,
based on a second, possible incomplete, reconstruction in which those dendrites
have been removed.

The "unninervated" version of the file does not have to include the soma or
axon.


This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021, 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""


import re
from pathlib import Path

import numpy as np
import pylibrary.tools.cprint as CP
from neuron import h

cprint = CP.cprint

file_full = Path(
    "/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c09/Morphology/VCN_c09_Full_MeshInflate.hoc"
)
print(" Full File found: ", file_full.is_file())
file_trim = Path(
    "/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c09/Morphology/VCN_c09_NoUninnervated.hoc"
)
print(" Trimmed File found: ", file_trim.is_file())
file_trim2 = Path(
    "/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c09/Morphology/VCN_c09_NoUninnervated_MeshInflate.hoc"
)

re_pts = re.compile(
    "\s*(pt3dadd\()([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)"
)
re_section = re.compile("\s*(sections\[)([0-9]*)\]\s*{")
re_access = re.compile("\s*(access)\s*(sections\[)([0-9]*)\]\s*")
re_append = re.compile("\s*(\w+)(.append\(\))")
re_connect = re.compile(
    "\s*(connect)\s*(sections\[)([0-9]*)\](\([0-9]*\)),\s*(sections\[)([0-9]*)\](\([0-9]*\))"
)
re_id = re.compile(
    "\s*(pt3dadd\()([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\) // SWC ID = ([0-9]*)"
)
re_objref = re.compile(r"\s*(objref) (\w+)\s*")


alwayskeeplist = [
    "soma",
    "Axon_Hillock",
    "Myelinated_Axon",
    "Unmyelinated_Axon",
    "Axon_Initial_Segment",
    "Axon_Heminode",
    "Axon_Node",
]

"""
keepids is a short list of the IDs we need to keep, and marked by a comment
with the section number
"""
keepids = [
    2003,
    2004,  # 106
    2010,
    2011,
    2012,
    2013,  # 117
    2005,
    2006,
    2007,
    2008,
    2009,  # 107
]  # always keep these IDs (has to do with the way the no-unninervated was reconstructed)


def get_ids(filename):
    ids = []
    last_section = None
    with open(Path(filename), "r") as fh:
        xf = fh.readlines()
    for line in xf:
        # n = re_access.match(line)  # the sections will not match...
        # if n is not None:
        #     last_section = int(n.groups()[2])
        m = re_id.match(line)
        if m is not None:
            ID = m.groups()[5]
            ids.append(int(ID))
    return ids
    # print(".. ", n)


def diff_files(file_full, file_trim, flag: str = "comment"):
    """
    compare annotated two hoc files and get the difference in the SWC IDs

    flag: must be "comment" or "delete"
        determines whether the output has the removed sections just commented
        out with // or actually deleted (not written)

    """
    assert flag in ["comment", "delete"]
    id1 = set(get_ids(file_full))
    id2 = set(get_ids(file_trim))

    nomatch = list(id1.difference(id2))
    # print("nomatch: ", sorted(nomatch))
    nd = np.diff(sorted(nomatch))
    print(
        f"# of entries: fullfile: {len(id1):d}, filetrim={len(id2):d}, nomatch={len(nomatch):d}"
    )
    lineout = []
    with open(Path(file_full)) as fh:
        xf = fh.readlines()
        sectionlines = []
        objrefs = []  # names of all objects found
        in_section = False
        appended_types = []
        removed_sections = []
        n3d_pts = 0
        remove_section = False
        id_list = []
        current_section = None
        for line in xf:
            if line.startswith("//"):
                lineout.append(line)  # echo comment lines
                continue
            if line.startswith("create sections"):
                lineout.append(line)
                continue
            sec = re_access.match(line)  # look for access statements
            if sec is not None:
                in_section = True  # start of section definition
                sectionlines.append(line)
                current_section = int(sec.groups()[2])
                continue

            if in_section:
                if line.startswith("sections["):
                    sectionlines.append(line)
                    continue
                ma = re_append.match(line)
                if ma is not None:
                    sectype = ma.groups()[0]
                    if sectype not in appended_types:  # keep track of types in use
                        appended_types.append(sectype)
                    sectionlines.append(line)
                    continue
                if re_connect.match(line) is not None:  # always add the connections
                    sectionlines.append(line)
                    continue
                m1 = re_id.match(line)
                if m1 is not None:  # line has an swc ID
                    if sectype in alwayskeeplist:
                        sectionlines.append(line)
                        n3d_pts += 1
                        continue
                    ID = int(m1.groups()[5])  # check the id
                    if ID in nomatch and ID not in keepids:
                        remove_section = True
                        id_list.append(ID)
                    else:
                        n3d_pts += 1
                    sectionlines.append(line)
                    continue

                if line.startswith("}"):  # closure
                    if remove_section:  # section flagged for removal
                        print(
                            f"Section[{current_section:4d}] of type {sectype:s} was removed, found in 'nomatch'"
                        )
                        print("    IDS were: ", id_list)
                        if flag == "comment":
                            sectionlines.append(line)
                            for i, sl in enumerate(sectionlines):
                                sectionlines[
                                    i
                                ] = f"//{sl:s}"  # comment them out. Could also just not add to lineout
                            lineout.append(sectionlines)
                        n3d_pts = 0
                        id_list = []
                        sectionlines = []  # clear the list of lines here
                        current_section = None
                        in_section = False  # end of section
                        remove_section = False
                    else:
                        sectionlines.append(line)
                        lineout.append(sectionlines)
                        n3d_pts = 0
                        sectionlines = []
                        in_section = False

            else:
                mobj = re_objref.match(line)
                if mobj is not None:
                    objrefs.append(mobj.groups()[1])
                lineout.append(line)

    with open(Path(file_trim2), "w") as fh:
        for l in lineout:
            if isinstance(l, list):
                for lx in l:
                    fh.write(lx)
            else:
                fh.write(l)
    return file_trim2


def delete_unparented_sections(file_trim2):
    ft = h.load_file(str(file_trim2))
    for i, section in enumerate(h.allsec()):
        sref = h.SectionRef(sec=section)
        psec = h.parent_section(0, sec=section)
        if not sref.has_parent():
            print("i: ", i, "  sec: ", section.name(), "has no parent")
            if i > 0:
                h.delete_section(sec=section)
        # if psec == 0.0:
        #     print(section.name())
        #     h.delete_section(sec=section)

    h.topology()  # to confirm that they do not exist...


if __name__ == "__main__":
    file_trim2 = diff_files(file_full, file_trim, flag="comment")
    delete_unparented_sections(file_trim2)
