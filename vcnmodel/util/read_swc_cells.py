"""
Read the swc soma reconstructions and place them into a 3-D space

"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as mpl
import swc_to_hoc

swcPath = Path('/Users/pbmanis/Desktop/Python', 'VCNModel', 'ASA', 'CellBodySWCs')

swcFiles = swcPath.glob('*.swc')

#print([f for f in swcFiles])

for f in swcFiles:
    SWC = swc_to_hoc.SWC(filename=f)
    print('topology for: ', f.name)
    SWC.write_hoc(Path(f.with_suffix('.hoc')))
     
    
    