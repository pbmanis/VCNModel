"""
List the active morphology files in our data set

Useful for checking which files are available and their current status.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021, 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

from datetime import datetime
from pathlib import Path

basepath = Path("VCN_Cells")

for cellnum in [2, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 24, 27, 29]:
    cellname = f"VCN_c{cellnum:02d}"
    cellpath = Path(basepath, cellname, "Morphology")
    mfiles = cellpath.glob("*.hoc")
    print("\n", cellname)
    for m in list(mfiles):
        dtime = datetime.fromtimestamp(m.stat().st_mtime)
        print(f"{str(m.name):32s}  {str(dtime):s}")
    if len(list(mfiles)) == 0:
        print(f"{' ':32s}  <No hoc files>")
