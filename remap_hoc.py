#!/usr/bin/env python
# encoding: utf-8
"""
remap_hoc.py

Created by Paul Manis on 2013-01-04.
Copyright (c) 2013 Paul B. Manis, Ph.D.. All rights reserved.

This function remaps a hoc file from 'axon[x]' to whatever is in the sectionLists.

It turns sections only found in the "section list" and replaces axon[k] with an appropriate 'swelling[n]'
A
ll is hardcoded...
"""

import sys
import os
import re
import numpy as np

#scaleFactor = 0.0384025 # convert pixels to microns for one file.
scaleFactor = [1.0, 1.0, 1.0, 1.0]
scaleFactor = [0.1*s for s in scaleFactor]
translateFlag = True

structures = ['soma', 'dend'] # ['tip', 'neck', 'swelling', 'branch', 'heminode', 'stalk', 'parentaxon']

def main():
    global xref
    # infile = 'MorphologyFiles/Calyx-S53Acvt3.hoc'
    # outfile = 'MorphologyFiles/Calyx-S53Acvt4.hoc'
    infile = 'MorphologyFiles/mainDenHOC.hoc'
    outfile = 'MorphologyFiles/mainDenHOC_cleaned.hoc'
   # infile = 'LC_neuromantic_scaled.hoc'
   # outfile = 'LC_nmscaled_cleaned.hoc'
    axonfind = re.compile('\{(?P<source>axon\[\d+\]) connect (?P<target>axon\[\d+\])\(0\), 1\}')
    connectfind = re.compile('\{axon\[(?P<source>\d+)\] connect axon\[(?P<target>\d+)\]\(0\), 1\}')
    accessfind = re.compile('\{access (?P<source>axon\[\d+\])\}(?P<comment>.*)')
    pt3dclearfind = re.compile('\{pt3dclear\(\)\}(?P<comment>.*)')
    pt3daddfind = re.compile('\{(?P<point>pt3dadd\((([-+]?\d*[\.]?\d*|\d+)[,]?[\s]*){1,4}\)){1,1}\}(?P<comment>.*)')
    points = re.compile('\s?\pt3dadd\((?P<point>((([-+]?\d*[\.]?\d*|\d+)[,]?[\s]*){1,4}))\)(?P<comment>.*)')
    anyaxon = re.compile('(?P<type>axon)\[(?P<number>\d+)\]')
    sre = {}
    smatch = {}
    xref = {}
    for s in structures:
        xref[s] = [] # will fill with list of axons, in order...
        sre[s] = re.compile('(?P<type>%s) = SectionList\(\)' % s) # top line of the group
        smatch[s] = re.compile('axon\[(?P<axon>\d+)\] %s\.append\(\)' % s)


    firstPoint = True
    # pt3dFlag = False
    inf = open(infile, 'r')
    zeropos = [0., 0., 0.]
    for lct, line in enumerate(iter(inf)):
        #print line[0:9]
        for s in structures:
            l = smatch[s].match(line)
            if l is not None:
                xref[s].append((int(l.group('axon')), len(xref[s]))) # build the substituion list
        mo  = points.match(line)
        if lct < 40:
            print 'mo, line: ', mo, line,
        if mo is not None:
            #mo = re.search(r"pt3dadd", line)
            nl = np.array(eval(mo.group('point')))
          #  nl = np.array(eval(line[8:]))
            if translateFlag:
                if firstPoint:
                    zeropos = nl[0:3].copy() # save the very first position
                    print 'translate to: ', zeropos
                    firstPoint = False
                nl[0:3] = nl[0:3] - zeropos
            nls = (nl * scaleFactor)
            nls = (nls)
            newline = '\tpt3dadd('
            for i, f in enumerate(nls):
                newline += str(f)
                if i < nls.shape[0]-1:
                    newline  += ', '
            newline += ')'+  '\n'
          #  print 'newline: ', newline,


    inf.close()
    # second pass: replace using the information in xref...
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    for cs in xref.keys(): # rewrite creation
        outf.write('create %s[%d]\n' % (cs, len(xref[cs])))
    for line in iter(inf):
        #print line[0:9]
        newline = line
        l = connectfind.match(line) # deal with the connections
        if l is not None:
            (ntype, n) = findmap(xref, int(l.group('source')))
            (mtype, m) = findmap(xref, int(l.group('target')))
            if ntype is not None and mtype is not None:
                newline = '{%s[%d] connect %s[%d](0), 1}\n' % (ntype, n, mtype, m)
                outf.write(newline)
                continue
        l = anyaxon.findall(line)
        if len(l) > 0:
            for p in l: # look for all the replacements.
                if submap(int(p[1])) == None:
                    print 'submap none:'
                    print newline
                    continue
                newline = anyaxon.sub(submap(int(p[1])), newline)
        if len(newline) == 1 and pt3dFlag:
            line = '}\n' # provide a closing brace
            pt3dFlag = False
        l = axonfind.match(newline)
        if l is not None:
           newline = '\t%s connect %s(1), 0\n' % (l.group('source'), l.group('target'))
           #line = '\tconnect %s(1), %s(0)\n' % (l.group('source'), l.group('target'))
        l = accessfind.match(newline)
        if l is not None:
            newline = '%s {\n' % (l.group('source')) #, l.group('comment'))
        l = pt3dclearfind.match(newline)
        if l is not None:
            newline = '\tpt3dclear()\n' # % (l.group('comment'))
        l = pt3daddfind.match(newline)
        if l is not None:
            newline = '\t' + l.group('point') + '\n'

#            line = '\t' + l.group('point') + l.group('comment') + '\n'
        #    print 'pt3daddfind done...'
            pt3dFlag = True # flag that we set a 3d point on this line
        else:
            pass
        #    print 'line failed to match pt3addfind: ', line

        # now scale if it's a pt3dadd.
        mo  = points.match(newline)
        if mo is not None:
            #mo = re.search(r"pt3dadd", line)
            nl = np.array(eval(mo.group('point')))
          #  nl = np.array(eval(line[8:]))
            if translateFlag:
                nl[0:3] = nl[0:3] - zeropos
            nls = (nl * scaleFactor)
            nls = (nls)
            newline = '\tpt3dadd('
            for i, f in enumerate(nls):
                newline += str(f)
                if i < nls.shape[0]-1:
                    newline  += ', '

#            line += ')'+ mo.group('comment')  + '\n' # may have a comment to follow
            newline += ')'# +  '\n' # may have a comment to follow
            outf.write(newline+'\n')
        else:
            pt3dFlag = False # no such flag.
            outf.write(newline)

    inf.close()
    outf.close()



def findmap(map, axno):
    for s in map:
        for (n,m) in map[s]:
            if axno == n:
                return(s, m) # return the new mapping.
    return(None, None)

def submap(axno):
    global xref
    for s in xref:
        for (n,m) in xref[s]:
            if int(axno) == n:
                return('%s[%d]' %(s, m) ) # return the new mapping.
    return(None)


if __name__ == '__main__':
    main()

