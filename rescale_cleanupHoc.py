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
scaleFactor = [1.0, 1.0, 1.0, 1.]


def main():
    infile = 'Calyx-68cvt2.hoc'
    outfile = 'Calyx-68cvt2_neurovisio_scaled.hoc'
    axonfind = re.compile('\{(?P<source>axon\[\d+\]) connect (?P<target>axon\[\d+\])\(0\), 1\}')
    accessfind = re.compile('\{access (?P<source>axon\[\d+\])\}(?P<comment>.*)')
    pt3dclearfind = re.compile('\{pt3dclear\(\)\}(?P<comment>.*)')
    pt3daddfind = re.compile('\{(?P<point>pt3dadd\((([-+]?\d*[\.]?\d+|\d+)[,]?){1,4}\)){1,1}\}(?P<comment>.*)')
    points = re.compile('\s?\pt3dadd\((?P<point>((([-+]?\d*[\.]?\d+|\d+)[,]?){1,4}))\)(?P<comment>.*)')

    pt3dFlag = False
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    for line in iter(inf):
        #print line[0:9]
        # do line replacements first
        #print len(line)
        if len(line) == 1 and pt3dFlag:
            line = '}\n' # provide a closing brace
            pt3dFlag = False
        l = axonfind.match(line)
        if l is not None:
           line = '\t%s connect %s(1), 0\n' % (l.group('source'), l.group('target'))
           #line = '\tconnect %s(1), %s(0)\n' % (l.group('source'), l.group('target'))
        l = accessfind.match(line)
        if l is not None:
            line = '%s {\n' % (l.group('source')) #, l.group('comment'))
        l = pt3dclearfind.match(line)
        if l is not None:
            line = '\tpt3dclear()\n' # % (l.group('comment'))
        l = pt3daddfind.match(line)
        if l is not None:
            line = '\t' + l.group('point') + '\n'
            
#            line = '\t' + l.group('point') + l.group('comment') + '\n'
        #    print 'pt3daddfind done...'
            pt3dFlag = True # flag that we set a 3d point on this line
        else:
            pass
        #    print 'line failed to match pt3addfind: ', line
        
        # now scale if it's a pt3dadd. 
        # print 'cleaned line: ', line
        mo  = points.match(line)
        if mo is not None:
        #    print 'mo: ', mo
        #    print mo.group()
        #    print 'mo group: ', mo.group('point')
        #    print 'new array: ', np.array(eval(mo.group('point')))
            #mo = re.search(r"pt3dadd", line)
            nl = np.array(eval(mo.group('point')))
          #  print 'original: ', line[8:]
          #  nl = np.array(eval(line[8:]))
            nls = (nl * scaleFactor)
          #  print 'scaled:    ', nls
            nls = (nls)
            line = '\tpt3dadd('
            for i, f in enumerate(nls):
                line += str(f)
                if i < nls.shape[0]-1:
                    line  += ', '
                
#            line += ')'+ mo.group('comment')  + '\n' # may have a comment to follow
            line += ')'+  '\n' # may have a comment to follow
         #   outf.write(newLine+'\n')
        else:
            pt3dFlag = False # no such flag.
        outf.write(line)
        #print line,
    inf.close()
    outf.close()


if __name__ == '__main__':
    main()

