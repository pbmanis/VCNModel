#!/usr/bin/env python
# encoding: utf-8
"""
rescaleHoc.py

Created by Paul Manis on 2013-01-04.
Copyright (c) 2013 Paul B. Manis, Ph.D.. All rights reserved.

This function rescales a hoc file. 
All is hardcoded... 
"""

import sys
import os
import re
import numpy as np

#scaleFactor = 0.0384025 # convert pixels to microns for one file.
scaleFactor = 1.0
def main():
    infile = 'Calyx-S53A_neurovisio.hoc'
    outfile = 'Calyx-S53A_neurovisio_scaled.hoc'
    
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    for line in iter(inf):
        #print line[0:9]
        mo = re.search(r"pt3dadd", line)
        if mo:
            print 'original: ', line[8:]
            nl = np.array(eval(line[8:]))
            nls = (nl * scaleFactor)
            print 'scaled:    ', nls
            nls = (nls)
            newLine = line[0:8] + '('
            for i, f in enumerate(nls):
                newLine += str(f)
                if i < nls.shape[0]-1:
                    newLine  += ', '
                
            newLine += ')'
            print 'new line:   ', newLine
            outf.write(newLine+'\n')
        else:
            outf.write(line)
    inf.close()
    outf.close()


if __name__ == '__main__':
    main()

