from __future__ import print_function
__author__ = 'pbmanis'

import sys
import os
from collections import OrderedDict
import pyqtgraph as pg
from cnmodel import cells
from cnmodel.decorator import Decorator
#from neuronvis import HocViewer

allcells = ['09', '11', '14', '17', '18', '19', '20', '21', '22']

def area(fn):
    filename = os.path.join('VCN_Cells', fn, 'Morphology', fn+'.hoc')
    post_cell = cells.Bushy.create(morphology=filename, decorator=Decorator,
            species='mouse',
            modelType='II', modelName='XM13')
    post_cell.distances()
    post_cell.computeAreas()
    secareas = {}
    #print 'areamap: ', post_cell.areaMap
    for am in post_cell.areaMap.keys():
        sectype = post_cell.get_section_type(am)
        if sectype is None:
            continue
        if sectype not in secareas.keys():
            secareas[sectype] = post_cell.areaMap[am]
        else:
            secareas[sectype] = secareas[sectype] + post_cell.areaMap[am]
    return secareas
    
    
if __name__ == '__main__':
    
    ar = OrderedDict()
    for fn in allcells:
        filename = 'VCN_c'+fn # os.path.join('VCN_Cells', 'VCN_c'+fn, 'Morphology', 'VCN_c'+fn+'.hoc')
        ar[fn] = area(filename)
    print('')
    hdrstr = ['Cell', 'Somatic area', 'Dendritic area', 'Ratio', 'Hillock Area', 'Unmyel Area', 'Myelin Area']
    hdrkeys = ['', 'soma', 'dendrite', 'ratio', 'hillock', 'unmyelinatedaxon', 'myelinatedaxon']
    dec = [0, 2, 2, 3, 2, 2, 2]  # define decimals for each column
    
    headers = ''.join('{:^12s}  '.format(hs) for hs in hdrstr)
    print( headers)
    for i, fn in enumerate(allcells):
        ar[fn]['ratio'] = ar[fn]['dendrite']/ar[fn]['soma']
        txt = ''
        for ik, k in enumerate(hdrkeys):
            if k == '':
                txt += '{:^12s}  '.format('VCN_c{0:2s}'.format(fn))
            else:
                txt += '{:>12.{:d}f}  '.format(ar[fn][k], dec[ik])
        print( txt)
        # print ('VCN_c{0:2s}:  Somatic area: {1:>8.2f}   Dendritic area: {2:>8.2f}  Ratio: {3:>5.3f}'
        #     .format(fn, ar[fn]['soma'], ar[fn]['dendrite'], ar[fn]['dendrite']/ar[fn]['soma']))
