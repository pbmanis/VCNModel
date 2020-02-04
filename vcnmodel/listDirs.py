from __future__ import print_function
__author__ = 'pbmanis'

import pickle
import pprint
import os
import sys

def listDirs():
    rootDir = 'Simulations/'
    for dirName, subdirList, fileList in os.walk(rootDir):
        print('Found directory: %s' % dirName)
        for fname in fileList:
            with open(dirName+fname, "r") as f:
                d = pickle.load(f)
                print('File: %s\n', fname)
                pprint.pprint(d['runInfo'], indent=4)


if __name__ == "__main__":
    listDirs()