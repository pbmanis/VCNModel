from pathlib import Path
import re
import numpy as np
from neuron import h
import pylibrary.tools.cprint as CP

cprint = CP.cprint

file_full = Path("/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c09/Morphology/VCN_c09_Full_MeshInflate.hoc")
print(" Full File found: ", file_full.is_file())
file_trim = Path("/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c09/Morphology/VCN_c09_NoUninnervated.hoc")
print(" Trimmed File found: ", file_trim.is_file())
file_trim2 = Path("/Volumes/Pegasus_002/VCN-SBEM-Data/VCN_Cells/VCN_c09/Morphology/VCN_c09_NoUninnervated_MeshInflate.hoc")

re_pts = re.compile('\s*(pt3dadd\()([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)')
re_section = re.compile('\s*(sections\[)([0-9]*)\]\s*{')
re_access = re.compile('\s*(access)\s*(sections\[)([0-9]*)\]\s*')
re_append = re.compile('\s*(\w+)(.append\(\))')
re_connect = re.compile('\s*(connect)\s*(sections\[)([0-9]*)\](\([0-9]*\)),\s*(sections\[)([0-9]*)\](\([0-9]*\))')
re_id = re.compile('\s*(pt3dadd\()([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\) // SWC ID = ([0-9]*)')
re_objref = re.compile(r'\s*(objref) (\w+)\s*')


alwayskeeplist =  ['soma', 'Axon_Hillock', 'Myelinated_Axon', 'Unmyelinated_Axon', 
                    'Axon_Initial_Segment', 'Axon_Heminode', 'Axon_Node']

"""
keepids is a short list of the IDs we need to keep, and marked by a comment
with the section number
"""
keepids = [2003, 2004,  # 106
           2010, 2011, 2012, 2013, #117
           2005, 2006, 2007, 2008, 2009,  # 107
       ]  # always keep these IDs (has to do with the way the no-unninervated was reconstructed)

def get_ids(filename):
    ids = []
    last_section = None
    with open(Path(filename), 'r') as fh:
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

id1 = set(get_ids(file_full))
id2 = set(get_ids(file_trim))

nomatch = list(id1.difference(id2))
print(sorted(nomatch))
nd = np.diff(sorted(nomatch))
print(len(id1), len(id2), len(nomatch))
lineout = []
with open(Path(file_full)) as fh:
    xf = fh.readlines()
    sectionlines = []
    objrefs = []  # names of all objects found
    in_section = False
    appended_types = []
    n3d_pts = 0
    p_nomatch = False
    for line in xf:
        if line.startswith("//"):
            lineout.append(line)  # echo comment lines
            continue
        if line.startswith("create sections"):
            lineout.append(line)
            continue
        sec = re_access.match(line)  # look for access statements
        if sec is not None:
            in_section = True   # start of section definition
            sectionlines.append(line)
            continue

        if in_section:
            if line.startswith("sections["):
                sectionlines.append(line)
            elif re_append.match(line) is not None:
                ma = re_append.match(line)
                sectype = ma.groups()[0]
                if sectype not in appended_types:
                    appended_types.append(sectype)
                sectionlines.append(line)
            elif re_connect.match(line) is not None:
                sectionlines.append(line)
            elif re_id.match(line) is not None:
                m1 = re_id.match(line)
                if m1 is None or sectype in alwayskeeplist:
                    sectionlines.append(line)
                    n3d_pts += 1
                else:
                    ID = int(m1.groups()[5])  # check the id
                    # if ID == 2014:  # a point in the soma
                    #     cprint("r", f"ID: {ID:5d}, {str(line):s}")
                    #     cprint("r", f" {str(ID in nomatch):s}, {sectype:s}")
                    if ID in nomatch and ID not in keepids:
                        # if not p_nomatch:
                        #     cprint("y", f"Found ID {ID:d} in nomatch: {str(nomatch):s}")
                        #     p_nomatch=True
                        # else:
                        #     cprint("y", f"Found ID {ID:d} in nomatch")
                        pass
                    else:
                        # if sectype == "soma":
                        #      cprint("m", "ID was ok to keep")
                        # else:
                        #     cprint("c", f"ID not soma, but kept, type = <{sectype:s}>")
                        n3d_pts += 1
                        sectionlines.append(line)
            elif line.startswith("}"):  # closure
                if n3d_pts == 0:  # no points
                    in_section = False  # end of section
                    print(f"Section of type = {sectype:s} was deleted as no pt3d found:")
                    print(sectionlines)
                    print()
                    sectionlines = [] # clear the list of lines here 
                    
                else:
                    sectionlines.append(line)
                    lineout.append(sectionlines)
                    n3d_pts = 0
                    sectionlines = []
                    in_section = False
                    
            
        else:
            mo = re_objref.match(line)
            if mo is not None:
                objrefs.append(mo.groups()[1])
            lineout.append(line)
                    
         #
        # m = re_id.match(line)
        # if m is None:
        #     lineout.append(line)  # no id on the line, so echo back
        # else:
        #     ID = int(m.groups()[5])  # check the id
        #     if ID in nomatch:
        #         line = '//' + line  # comment the line out
        #         lineout.append(line)
        #     else:
        #         lineout.append(line)  # just print it
        #
        # if not in_section:
        #     lineout.append(line)

with open(Path(file_trim2), 'w') as fh:
    for l in lineout:
        if isinstance(l, list):
            for lx in l:
                fh.write(lx)
        else:
            fh.write(l)

ft = h.load_file(str(file_trim2))
h.topology()
# print(objrefs)
# print(appended_types)