"""
List the active morphology files in our data set

Useful for checking which files are available and their current status.
"""

from pathlib import Path
from datetime import datetime

basepath = Path('VCN_Cells')

for cellnum in [2, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 24, 27, 29]:
    cellname = f"VCN_c{cellnum:02d}"
    cellpath = Path(basepath, cellname, 'Morphology')
    mfiles = cellpath.glob('*.hoc')
    print('\n', cellname)
    for m in list(mfiles):
        dtime = datetime.fromtimestamp(m.stat().st_mtime)
        print(f"{str(m.name):32s}  {str(dtime):s}")
    if len(list(mfiles)) == 0:
        print(f"{' ':32s}  <No hoc files>")
