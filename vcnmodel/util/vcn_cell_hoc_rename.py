"""
copy the ORIGINAL VCN hoc files to a new name only if the new name does not already exist (prevent overwrite)

"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as mpl
import neuronvis.swc_to_hoc as swc_to_hoc
import toml
config = toml.load(open("wheres_my_data.toml", "r"))

def main():
    celln = [2, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 27, 29]
    cellPath = Path(config['codeDirectory'], 'VCN_Cells')

    VCNCells = list(cellPath.glob('VCN_c*'))

    print(VCNCells)
    if len(VCNCells) == 0:
        print('No files or directory no longer exists')
        exit()
    
    for f in VCNCells:
        hocPath = Path(cellPath, f.name, 'Morphology', f.name+'.hoc')
        if hocPath.exists():
            destination = Path(cellPath, f.name, 'Morphology', f.name+'_original.hoc')
            print('old: ', hocPath, '\nnew: ', destination)
            with destination.open(mode='xb') as fid:
                fid.write(hocPath.read_bytes())

     
if __name__ == '__main__':
     main()   
    