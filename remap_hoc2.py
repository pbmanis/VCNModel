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
import pprint
import neuronvis.hoc_reader as HR

#scaleFactor = 0.0384025 # convert pixels to microns for one file.
scaleFactor = np.array([1.0, 1.0, 1.0, 1.0]) # multiply incoming by this. isotropic.
#scaleFactor =  np.array([48., 48., 48., 2.0*48.])/1000. # for data from Morehead/Spirou

translateFlag = False

structures = ['soma', 'dend'] # ['tip', 'neck', 'swelling', 'branch', 'heminode', 'stalk', 'parentaxon']

cnxns = {0: None,
         1: [0, [0, 0, 7]],
         2: [0, [8, 8, -7]],
         3: [0, [-8, -8, -7]],
         4: [0, [0, 0, -17]],
         5: [4, [0, 0, 0]],
         6: [4, [0, 0, 0]],
         7: [6, [0, 0, 0]],
         8: [6, [0, 0, 0]],
        }


def main():
    global xref
#    infile = 'Calyx-68cvt2.hoc'
#    outfile = 'Calyx-68cvt4.hoc'
    path = 'MorphologyFiles/'
   # infile = 'LC_neuromantic_scaled.hoc'
   # outfile = 'LC_nmscaled_cleaned.hoc'
    #infile='wholeThing.hoc'
    #outfile='wholeThing_cleaned.hoc'
    infile = 'MNTB_Cell2_cleaned_V2.hoc' # _cleaned_V2.hoc'
    outfile = 'MNTB_Cell2_cleaned_V2a.hoc'
    infile = os.path.join(path, infile)
    outfile = os.path.join(path, outfile)
    hrf = HR.HocReader(infile)
    names = {0: 'soma[0]'}
    for i in range(0, 8):

        names[i+1] = 'dend[%d]' % i
   # print names
    geom =  hrf.get_geometry()
    clist = geom[1]
    plist = geom[0]
    rlist = np.zeros(len(plist))
    dlist = {}
    zero = True
#    print 'plist[0][0]', plist[0][0]
    for x in range(len(plist)):
        if zero:
            p0 = plist[x][0]
            zero = False
        plist[x][0] = plist[x][0] - p0
        plist[x][0] = plist[x][0]/scaleFactor[0:3]
        rlist[x] = plist[x][1] # radius list
        print plist[x][2], x # element list
        dlist[plist[x][2]].append(x)]
        print dlist

    exit()

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
    structs = {}
    for s in structures:
        xref[s] = [] # will fill with list of axons, in order...
        sre[s] = re.compile('(?P<type>%s) = SectionList\(\)' % s) # top line of the group
        smatch[s] = re.compile('axon\[(?P<axon>\d+)\] %s\.append\(\)' % s)
        structs[s] = re.compile('%s\[(?P<number>\d+)\]*' % s)



    dendmap = {'dend[0]': (0., 0., -17.56)}
    transmap = []
    firstPoint = True
    # pt3dFlag = False
    print 'Processing: ', infile
    inf = open(infile, 'r')
    zeropos = [0., 0., 0.]
    allstructs = {}

    cpoint = [0., 0., 0.] # current point
    laststruct = None
    thisstruct = None
    for lct, line in enumerate(iter(inf)):
        #print line[0:9]
        for s in structures:
            l = smatch[s].match(line)
            if l is not None:
                xref[s].append((int(l.group('axon')), len(xref[s]))) # build the substituion list
            xl = structs[s].match(line)
            if xl is not None:
               # print 'this struct: ', s
                n = int(xl.group('number'))
                thisstruct = "%s[%d]" % (s, n)
                if thisstruct != laststruct:
                   laststruct = thisstruct

                if thisstruct not in allstructs.keys():
                    print 'new: %s' % thisstruct
                    allstructs[thisstruct] = [None, None]
                    allstructs[thisstruct][1] = cpoint


        mo  = points.match(line)
        #if lct < 40:
        #    print 'mo, line: ', mo, line,
        if mo is not None:
            #mo = re.search(r"pt3dadd", line)
            nl = np.array(eval(mo.group('point')))
            if translateFlag:
                if firstPoint:
                    zeropos = nl[0:3].copy() # save the very first position
                    print 'translate to: ', zeropos
                    firstPoint = False
                nl[0:3] = nl[0:3] - zeropos
            nls = (nl * scaleFactor)
            if allstructs[thisstruct][0] is None:
                allstructs[thisstruct][0] = nl[0:3] # get starting position of each process
            cpoint = nl[0:3]
            nls = (nls)
            newline = '\tpt3dadd('
            for i, f in enumerate(nls):
                newline += str(f)
                if i < nls.shape[0]-1:
                    newline  += ', '
            newline += ')'+  '\n'
          #  print 'newline: ', newline,

    pprint.pprint(allstructs)
    inf.close()
    # second pass: replace using the information in xref...
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
#    for cs in xref.keys(): # rewrite creation
#        outf.write('create %s[%d]\n' % (cs, len(xref[cs])))
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

