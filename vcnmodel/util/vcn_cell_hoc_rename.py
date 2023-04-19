"""
utility to copy the ORIGINAL VCN hoc files to a new name only if the new name
does not already exist (prevent overwrite)

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""

from pathlib import Path

from vcnmodel.util.get_data_paths import get_data_paths



def main():
    config = get_data_paths()

    celln = [2, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 27, 29]
    cellPath = Path(config["disk"], config["baseDataDirectory"], config["cellDataDirectory"])

    VCNCells = list(cellPath.glob("VCN_c*"))

    print(VCNCells)
    if len(VCNCells) == 0:
        print("No files or directory no longer exists")
        exit()

    for f in VCNCells:
        hocPath = Path(cellPath, f.name, "Morphology", f.name + ".hoc")
        if hocPath.exists():
            destination = Path(cellPath, f.name, "Morphology", f.name + "_original.hoc")
            print("old: ", hocPath, "\nnew: ", destination)
            # with destination.open(mode="xb") as fid:
            #     fid.write(hocPath.read_bytes())


if __name__ == "__main__":
    main()
