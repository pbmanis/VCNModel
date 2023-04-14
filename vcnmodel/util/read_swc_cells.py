"""
Read the swc soma reconstructions and genereate hoc files.

This data is not in the shared database.

"""

from pathlib import Path
import vcnmodel.util.swc_to_hoc as swc_to_hoc
from vcnmodel.util.get_data_paths import get_data_paths

config = get_data_paths()

swcPath = Path(config["disk"], config['baseMorphologyDirectory'], 'ASA', 'CellBodySWCs')

swcFiles = swcPath.glob('*.swc')

def main():

    for f in swcFiles:
        SWC = swc_to_hoc.SWC(filename=f)
        print('topology for: ', f.name)
        # SWC.write_hoc(Path(f.with_suffix('.hoc')))
     

if __name__ == "__main__":
    main()    
    