"""
Read a hoc file, and rescale x,y,z with an anisotropic set of factors
Designate an output file
"""
from __future__ import print_function
import sys
import os
import re
import numpy as np
import os.path

sf = 0.02
minimum_diameter = 0.25
scaleFactor = np.array([sf, sf, sf, sf])
path='VCN_Cells/VCN_c18/Morphology/'

def rescale_hoc_aniso(infile=None, outfile='temp.hoc'):
    assert infile is not None
#    path = 'MorphologyFiles/'
    infile = os.path.join(path, infile)
    outfile = os.path.join(path, outfile)

    points = re.compile('(?P<space>\s*)(?P<command>pt3dadd\()(?P<point>(([-+]?\d*[\.]?\d+|\d+)[,]?[\s*]?){1,4})\)(?P<comment>.*)')

    pt3dFlag = False
    inf = open(infile, 'r')
    print(infile)
    print(inf)
    outf = open(outfile, 'w')
    for line in iter(inf):
        
        # now scale if it's a pt3dadd. 
        mo  = points.match(line)
        if mo is not None:
            #mo = re.search(r"pt3dadd", line)
            nl = np.array(eval(mo.group('point')))
          #  nl = np.array(eval(line[8:]))
            print('nl : ', nl)
            nls = (nl * scaleFactor)
            nls = (nls)
            line = '\tpt3dadd('
            for i, f in enumerate(nls):
                if i == 3 and f < minimum_diameter:
                    f = minimum_diameter
                line += str(f)
                if i < nls.shape[0]-1:
                    line  += ', '
            print('new nls: ', nls)
#            line += ')'+ mo.group('comment')  + '\n' # may have a comment to follow
            line += ')'+  '\n' # may have a comment to follow
         #   outf.write(newLine+'\n')
        else:
            pt3dFlag = False # no such flag.
            #print 'did not match any points'
        outf.write(line)
        #print line,
    inf.close()
    outf.close()
    
if __name__ == '__main__':
    rescale_hoc_aniso(infile='gbc18_w_axon.hoc', outfile='gbc18_w_axon_rescaled.hoc')
    