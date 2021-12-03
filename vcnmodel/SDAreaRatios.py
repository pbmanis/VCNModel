#!/usr/bin/python
"""
Compute soma and dendrite areas and ratios"""
import sys
import numpy as np
from pathlib import Path
from collections import OrderedDict
import pyqtgraph as pg
from cnmodel import cells
from cnmodel.decorator import Decorator
import adjust_areas
import toml
config = toml.load(open("wheres_my_data.toml", "r"))

AdjA = adjust_areas.AdjustAreas()

gradeA = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]

all_cell_nos = range(1, 32)
allcells = [f"{c:02d}" for c in all_cell_nos]
# print(allcells)
__author__ = "pbmanis"

def sum_area(areamap):
    """
    areamap is post_cell.areaMap[sectiontype]
    is still a dict{ 'section[x]': floatingareavalue}
    """
    area = 0.0
    for x in areamap:
        area += float(areamap[x])
    return area


def one_sectype_area(amap):
    """
    Return the total area of one section type by name
    """
    area = 0.0
    area = np.sum([amap[s] for s in amap])
    # print(area)
    return area


def area(fn):
    """
    Get the area information for a cell
    specified as a hoc file
    """
    filename = Path(
        config['cellDataDirectory'],
        fn,
        "Morphology",
        fn + "_Full.hoc",
    )
    post_cell = cells.Bushy.create(
        morphology=str(filename),
        # decorator=Decorator,
        species="mouse",
        modelType="II",
        modelName="XM13",
    )
    # post_cell.list_sections()
    post_cell.distances()
    post_cell.computeAreas(source='seg')
    secareas = {}
    # print('soma: ', post_cell.areaMap['soma'])
    for am in list(post_cell.areaMap.keys()):
        # print('am: ', am)
        if am not in list(secareas.keys()):
            secareas[am] = one_sectype_area(post_cell.areaMap[am])
        else:
            secareas[am] += one_sectype_area(post_cell.areaMap[am])
    # for s in secareas.keys():
    #     print(f"    first pass {s:>24s} = {secareas[s]:.3f}")

    # now collapse dendrites...
    secareas["somabysegment"] = 0.0
    secareas["dendrite"] = 0.0
    # secareas['soma'] = 0.
    for am in list(post_cell.areaMap.keys()):
        # if am == "soma":
        #     secareas["soma"] += one_sectype_area(post_cell.areaMap[am])
        # print('am: ', am)
        if am.startswith("Dendr") or am.endswith("Dendrite") or am.startswith("dend"):
            # print("adding dendr area: ", am, one_sectype_area(post_cell.areaMap[am]))
            secareas["dendrite"] += one_sectype_area(post_cell.areaMap[am])
    # print(f"{filename.name:s}   secareas: {str(secareas):s}")
    return secareas


if __name__ == "__main__":

    ar = OrderedDict()
    for i, fn in enumerate(allcells):
        filename = (
            "VCN_c" + fn
        )  # os.path.join('VCN_Cells', 'VCN_c'+fn, 'Morphology', 'VCN_c'+fn+'.hoc')
        if i + 1 in gradeA:
            print("getting area for fn: ", fn)
            ar[fn] = area(filename)
        else:
            ar[fn] = {"soma": np.nan, "dendrite": np.nan}
    print("")
    hdrstr = [
        "Cell",
        "Somatic area",
        "Dendritic area",
        "Ratio",
        "Hillock Area",
        "Unmyel Area",
        "Myelin Area",
    ]
    hdrkeys = [
        "",
        "soma",
        "dendrite",
        "ratio",
        "hillock",
        "unmyelinatedaxon",
        "myelinatedaxon",
    ]
    dec = [0, 2, 2, 3, 2, 2, 2]  # define decimals for each column

    headers = "".join("{:^12s}  ".format(hs) for hs in hdrstr)
    # print(ar.keys())
    # print(headers)
    for i, fn in enumerate(ar.keys()):
        # print(ar[fn]["soma"], ar[fn]["dendrite"])
        # print(ar[fn].keys())
        ar[fn]["ratio"] = ar[fn]["dendrite"] / ar[fn]["soma"]
        txt = "{:^12s}  ".format("VCN_c{0:2s}".format(fn))
        if i + 1 == 29:
            print(f"{fn:^12s}.D02")
        for ik, k in enumerate(hdrkeys):
            if k not in list(ar[fn].keys()):
                continue
                # txt += "{:>12.{:d}f}  ".format(0.0, dec[ik])
            # elif k == "somabysegment":
            #     txt += "{:>12.{:d}f}  ".format(ar[fn][k], dec[ik])
            elif k == "":
                continue
            else:
                txt += f"{ar[fn][k]:>12.3f}  "  # .format(ar[fn][k], dec[ik])
        print(txt)
        if i + 1 == 24:
            print(f"{fn:^12s}.D01")
        # print ('VCN_c{0:2s}:  Somatic area: {1:>8.2f}   Dendritic area: {2:>8.2f}  Ratio: {3:>5.3f}'
        #     .format(fn, ar[fn]['soma'], ar[fn]['dendrite'], ar[fn]['dendrite']/ar[fn]['soma']))
