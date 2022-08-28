#!/usr/bin/env python3
'''
Created on Aug 9, 2017
@author: arnon

The MIT License (MIT)

Copyright (c) 2015 Acrisel

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Original: https://github.com/Acrisel/references/edit/master/osutils/ptree.py

Modified 8 Jul 2022 pbmanis
use pathlib, report file sizes

'''

import os  # need .walk, not provided by pathlib
from pathlib import Path

def realname(path, root=None):
    ''' joins root with path, if root is provided.
        Then check is it is a symlink.  If it is, return 
        a string representing the link.  Otherwise, return 
        basename or path.
    '''
    path=Path(path)
    if root is not None:
        path=Path(root, path)
    result=path.name
    if path.is_symlink():
        realpath=path.readlink()
        result= '%s -> %s' % (path.name, realpath)
    return result

def ptree(startpath, depth=-1): 
    ''' prints directory tree in 'tree' structure.
    
    Args:
        startpath: root path to start
        depth: depth of tree to print; default: -1 which signals not limit
    '''
    linesout = ""
    prefix=0
    startpath = str(startpath)
    if startpath != '/':
        if startpath.endswith('/'): startpath=startpath[:-1]
        prefix=len(startpath)  
    for root, dirs, files in os.walk(startpath):
        level = root[prefix:].count(os.sep)
        if depth >-1 and level > depth: continue
        indent=subindent =''
        if level > 0:
            indent = '|   ' * (level-1) + '|-- '
        subindent = '|   ' * (level) + '|-- '
        linesout += f"{indent:s}{realname(root)}/\n"
        # print dir only if symbolic link; otherwise, will be printed as root
        for d in dirs:
            if Path(root, d).is_symlink():
                linesout += f"{subindent:s}{realname(d, root=root):s}\n"
        for f in files:
            size = Path(root, f).stat().st_size
            if size > 1e3 and size < 1e6:
                s_size = f"{size/1e3:.2f} kB"
            elif size > 1e6 and size < 1e9:
                s_size = f"{size/1e6:.2f} MB"
            elif size > 1e9:
                s_size = f"{size/1e9:.1f} GB"
            else:
                s_size = f"{int(size):d} B"
            linesout += f"{subindent}{realname(f, root=root)}  [{s_size:s}]\n"
    return linesout

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="""prints directory tree""")
    parser.add_argument('--level', '-l', type=int, dest='depth', default=6, help='depth of tree to print')
    parser.add_argument('startpath', type=str, help='path to stating directory')
    args = parser.parse_args()
    argsd=vars(args)
    lo = ptree(**argsd)
    print(lo)