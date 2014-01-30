__author__ = 'pbmanis'

import matplotlib.pylab as MP
import pickle
import pylibrary.PlotHelpers as PH
import pprint
import os
import sys

def listDirs():
    rootDir = 'Canonical/'
    for dirName, subdirList, fileList in os.walk(rootDir):
        print('Found directory: %s' % dirName)
        for fname in fileList:
            with open(dirName+fname, "r") as f:
                d = pickle.load(f)
                print 'File: %s\n', fname
                pprint.pprint(d['runInfo'], indent=4)


if __name__ == "__main__":
    listDirs()