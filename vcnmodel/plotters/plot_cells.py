"""

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""


import matplotlib.pyplot as mpl
from pathlib import Path
import pylibrary.plotting.plothelpers as PH
import matplotlib.image as mpimg
import numpy as np
import toml



# In[63]:

def plot_cells():
    with open("wheres_my_data.toml", "r") as fh:
        config = toml.load(fh)
    
    fn = Path(config['baseDataDirectory'], 'VCN-CellImages-1-2020/Bushy_Cells')
    print(fn)
    fns = fn.glob('*.png')
    pngs = list(fns)
    npngs = len(pngs)
    print(npngs)
    rows, cols = PH.getLayoutDimensions(npngs)
    P = PH.regular_grid(rows, cols, order='rowsfirst', figsize=(10.0, 10.0), 
                        showgrid=False, verticalspacing=0.01, horizontalspacing=0.01, 
                        margins={'bottommargin': 0.1, 'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.03}, 
                        labelposition=(0.0, 0.0), parent_figure=None, panel_labels=None, )

    ext = 100
    axarray = P.axarr.ravel()
    for i, p in enumerate(pngs):
        img=mpimg.imread(str(p))
        axarray[i].imshow(img,  aspect='equal')
        axarray[i].set_title(str(p.name).replace('_', '\_'))
        ylim = axarray[i].get_ylim()
    #     print('xlim: ', ylim)
    # xlim for these ata are 0, 2500
        if 'c10' in str(p):
            axarray[i].set_xlim(200, 2000)
        else:
            axarray[i].set_xlim(800, 1700)
    #     axarray[i].set_ylim(1200, 300)
    
        PH.noaxes(axarray[i])
    for j in range(i, rows*cols):
        PH.noaxes(axarray[j])
    mpl.show()


if __name__ == '__main__':
    plot_cells()




