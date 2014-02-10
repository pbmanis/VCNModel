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

structures = ['tip', 'neck', 'swelling', 'branch', 'heminode', 'stalk', 'parentaxon']

def main():
    global xref
    # infile = 'MorphologyFiles/Calyx-S53Acvt3.hoc'
    # outfile = 'MorphologyFiles/Calyx-S53Acvt4.hoc'
    infile = 'MorphologyFiles/Calyx-68cvt2.hoc'
    outfile = 'MorphologyFiles/Calyx_68cvt4.hoc'
   # infile = 'LC_neuromantic_scaled.hoc'
   # outfile = 'LC_nmscaled_cleaned.hoc'
    axonfind = re.compile('\{(?P<source>axon\[\d+\]) connect (?P<target>axon\[\d+\])\(0\), 1\}')
    connectfind = re.compile('\{axon\[(?P<source>\d+)\] connect axon\[(?P<target>\d+)\]\(0\), 1\}')
    accessfind = re.compile('\{access (?P<source>axon\[\d+\])\}(?P<comment>.*)')
    pt3dclearfind = re.compile('\{pt3dclear\(\)\}(?P<comment>.*)')
    pt3daddfind = re.compile('\{(?P<point>pt3dadd\((([-+]?\d*[\.]?\d+|\d+)[,]?){1,4}\)){1,1}\}(?P<comment>.*)')
    points = re.compile('\s?\pt3dadd\((?P<point>((([-+]?\d*[\.]?\d+|\d+)[,]?){1,4}))\)(?P<comment>.*)')
    anyaxon = re.compile('(?P<type>axon)\[(?P<number>\d+)\]')
    sre = {}
    smatch = {}
    xref = {}
    for s in structures:
        xref[s] = [] # will fill with list of axons, in order...
        sre[s] = re.compile('(?P<type>%s) = SectionList\(\)' % s) # top line of the group
        smatch[s] = re.compile('axon\[(?P<axon>\d+)\] %s\.append\(\)' % s)


    # pt3dFlag = False
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    for line in iter(inf):
        #print line[0:9]
        for s in structures:
            l = smatch[s].match(line)
            if l is not None:
                xref[s].append((int(l.group('axon')), len(xref[s]))) # build the substituion list


    inf.close()
    # second pass: replace using the information in xref...
    inf = open(infile, 'r')
    # outf = open(outfile, 'w')
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

