from __future__ import print_function
__author__ = 'pbmanis'

import sys
from pathlib import Path
from collections import OrderedDict
import pyqtgraph as pg
from cnmodel import cells
from cnmodel.decorator import Decorator
#from neuronvis import HocViewer

allcells = ['08', '09', '11', '14', '16', '17', '18', '19', '20', '21', '22']


def sum_area(areamap):
    # areamap is post_cell.areaMap[sectiontype]
    # is still a dict{ 'section[x]': floatingareavalue}
    
    area = 0.
    for x in areamap:
        area += float(areamap[x])
    return area



def sum_area(areamap):
    # areamap is post_cell.areaMap[sectiontype]
    # is still a dict{ 'section[x]': floatingareavalue}
    
    area = 0.
    for x in areamap:
        area += float(areamap[x])
    return area


def area(fn):
    filename = Path('VCN_Cells', fn, 'Morphology', fn+'.hoc')
    post_cell = cells.Bushy.create(morphology=str(filename), decorator=Decorator,
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
            print('new type: %s  area: %f' % (sectype, post_cell.areaMap[am]))
        else:
            secareas[sectype] = secareas[sectype] + post_cell.areaMap[am]
            print('     type: %s  area: %f' % (sectype, secareas[sectype]))
    for am in post_cell.areaMap2.keys():
        sectype = post_cell.get_section_type(am)
        if sectype == 'soma':
            secareas['somabysegment'] = post_cell.areaMap2[am]

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
    for i, fn in enumerate(ar.keys()):

        ar[fn]['ratio'] = ar[fn]['dendrite']/ar[fn]['soma']
        txt = '{:^12s}  '.format('VCN_c{0:2s}'.format(fn))
        for ik, k in enumerate(hdrkeys):
            if k not in list(ar[fn].keys()):
                txt += '{:>12.{:d}f}  '.format(0., dec[ik])
            elif k == 'somabysegment':
                txt += '{:>12.{:d}f}  '.format(ar[fn][k], dec[ik])
            elif k == '':
                continue
            else:
                txt += '{:>12.{:d}f}  '.format(ar[fn][k], dec[ik])
        print( txt)
        # print ('VCN_c{0:2s}:  Somatic area: {1:>8.2f}   Dendritic area: {2:>8.2f}  Ratio: {3:>5.3f}'
        #     .format(fn, ar[fn]['soma'], ar[fn]['dendrite'], ar[fn]['dendrite']/ar[fn]['soma']))
