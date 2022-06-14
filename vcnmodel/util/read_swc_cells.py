"""
Read the swc soma reconstructions and genereate hoc files.

"""

from pathlib import Path
import vcnmodel.util.swc_to_hoc as swc_to_hoc
import toml
config = toml.load(open("wheres_my_data.toml", "r"))

swcPath = Path(config['baseMorphologyDirectory'], 'ASA', 'CellBodySWCs')

swcFiles = swcPath.glob('*.swc')

def main():

    for f in swcFiles:
        SWC = swc_to_hoc.SWC(filename=f)
        print('topology for: ', f.name)
        # SWC.write_hoc(Path(f.with_suffix('.hoc')))
     

if __name__ == "__main__":
    main()    
    