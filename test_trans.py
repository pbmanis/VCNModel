from parse_hoc import ParseHoc
import re
from pathlib import Path
pathn = Path('/Users/pbmanis/Desktop/Python/VCNModel/VCN_Cells/')
cellname = 'VCN_c08'
refpts = {'VCN_c11': 'first', 'VCN_c08': 'first', 'VCN_c09': 'last', 'VCN_c14': 'first', 'VCN_c16': 'first',
            'VCN_c17': 'first',  'VCN_c18': 'last', 'VCN_c19': 'first', 'VCN_c20': 'first',
            'VCN_c21': 'first', 'VCN_c22': 'first'}
hfname = Path(pathn, cellname, 'Morphology', cellname+'.hoc')
refname = Path(pathn, cellname, 'Morphology', cellname+'_original.hoc')

PH = ParseHoc(fn=hfname)
PH.print_sections(fref=refname, refpoint=refpts[cellname])
