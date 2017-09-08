__author__ = 'pbmanis'

import sys
import os
import pyqtgraph as pg
from cnmodel import cells
from cnmodel.decorator import Decorator
from neuronvis import HocViewer

def area(fn):
    filename = os.path.join('VCN_Cells', fn, 'Morphology', fn+'.hoc')
    post_cell = cells.Bushy.create(morphology=filename, decorator=Decorator,
            species='mouse',
            modelType='mGBC')
    post_cell.set_distances()
    post_cell.computeAreas()
    secareas = {}
    print 'areamap: ', post_cell.areaMap
    for am in post_cell.areaMap.keys():
        sectype = post_cell.get_section_type(am)
        if sectype is None:
            continue
        if sectype not in secareas.keys():
            secareas[sectype] = post_cell.areaMap[am]
        else:
            secareas[sectype] = secareas[sectype] + post_cell.areaMap[am]
    print secareas
    print 'Ratio dendrite/soma area: ', secareas['dendrite']/secareas['soma']
    
    
if __name__ == '__main__':
    
    for fn in ['09', '11', '14', '17', '18', '19', '20', '21', '22']:
        filename = os.path.join('VCN_Cells', fn, 'Morphology', fn+'.hoc')
        area(filename)
