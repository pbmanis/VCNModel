"""
Read the swc soma reconstructions and genereate hoc files.

"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as mpl
import swc_to_hoc
import toml
config = toml.load(open("wheres_my_data.toml", "r"))

swcPath = Path(config['baseMorphologyDirectory'], 'ASA', 'CellBodySWCs')

swcFiles = swcPath.glob('*.swc')


for f in swcFiles:
    SWC = swc_to_hoc.SWC(filename=f)
    print('topology for: ', f.name)
    # SWC.write_hoc(Path(f.with_suffix('.hoc')))
     
    
    