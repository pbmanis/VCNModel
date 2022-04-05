from __future__ import print_function

"""
Manage the vector strength summary files as .p files
Creates a dict with rows (cells) and columns (experimental manipulations)
Adds data to file, or prints results from the file

Usage:
import vspfile

routines:
vspfile.init_vsp()
vspfile.add_data()
vspfile.print_vsp()

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2014- Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
import os
import pickle
from collections import OrderedDict

import numpy as np

vspfilename = "VS.p"


def init_vsp(cells, runtypes):
    """
    init the vector strength summmary file
    If the file exists, we only add cells/runtypes to it
    If the file does not exist, we create an empty set of dicts.
    """
    if os.path.isfile(vspfilename):  # do not create if it already exists
        # check to see if all the cells and runtypes are in the table
        fh = open(vspfilename, "rb")  # set file handle to write  - open, and append
        vsp = pickle.load(fh)  # get data for the table
        fh.close()
        for c in cells:
            cn = "VCN_c%02d" % c
            if cn not in vsp.keys():  # add the cell and the runtype keys
                vsp[cn] = OrderedDict()
                for r in runtypes:
                    vsp[cn][r] = np.nan
            # check the runtypes for the existing cell
            else:
                for r in runtypes:
                    if r not in vsp[cn].keys():
                        vsp[cn][r] = np.nan

    else:  # create empty table with all cellkeys and runtype keys for each cell
        vsp = OrderedDict()  # create /initialize summary dict
        for c in cells:
            cn = "VCN_c%02d" % c
            vsp[cn] = {}
            for r in runtypes:
                vsp[cn][r] = np.nan

    fh = open(
        vspfilename, "wb"
    )  # set file handle to write  - open, which will make an empty file
    pickle.dump(vsp, fh)
    fh.close()


def add_data(rowname, colname, value):
    """
    Add data to the file. No check to see whether row'colum dict values already exist
    """
    fh = open(vspfilename, "rb")  # set file handle to write  - open, and append
    vsp = pickle.load(fh)  # get data for the table
    fh.close()
    vsp[rowname][colname] = value  # update entry for this cell in the table
    fh = open(vspfilename, "wb")  # rewrite
    pickle.dump(vsp, fh)  # save it.
    fh.close()


def print_vsp():
    """
    Print out the VS data file as a tab-formatted text that can be copied and imported
    into excel or prism, etc.
    """
    fh = open(vspfilename, "rb")
    vsp = pickle.load(fh)
    cells = vsp.keys()
    #    print cells
    runtypes = vsp[cells[0]].keys()
    #    print runtypes
    print("-" * 80)
    print(" " * 12, end="")
    for r in runtypes:
        print("{:>12s}\t".format(r), end="")
    print("")
    for c in vsp.keys():
        print("{0:<12s}\t".format(c), end="")
        for r in runtypes:
            if not np.isnan(vsp[c][r]):
                print("{0:>12.3f}\t".format(vsp[c][r]), end="")
            else:
                print("{0:>8s}\t".format("nan"), end="")
        print("")
    print("-" * 80)
    print("")


if __name__ == "__main__":
    """
    If called from command line, then just print the current file
    """
    print_vsp()
