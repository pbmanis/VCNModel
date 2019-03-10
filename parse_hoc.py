import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as mpl
import pandas as pd
import argparse
import re
import datetime

from neuron import h

re_pts = re.compile('\s*(pt3dadd\()([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)')
re_section = re.compile('\s*(sections\[)([0-9]*)\]\s*{')
re_append = re.compile('\s*([a-z]*)(.append\(\))')

"""
Parse a hoc file

"""
class ParseHoc(object):
    def __init__(self, fn=None):
        self.objrefs = {}
        self.fn = fn
        self.sections = []
    
    def valid_object(self, objref):
        if objref in ['soma', 'dendrite', 'hillock', 'unmyelinatedaxon', 'myelinatedaxon']:
            return True
        else:
            return False
            
    def read_hoc(self):
        lines = []
        with open(self.fn, 'r') as fh:
            for cnt, line in enumerate(fh):
                lines.append(line)         
        self.hoctxt = ' '.join(l for l in lines)
        print(self.hoctxt)
    
    def get_ref_point(self, fref, refpoint='first'):
        if fref is not None:
            self.fref = fref
        else:
            exit()
        xyzr = None
        base_xyz = None
        current_section_type = None
        got_first_point = False
        with open(self.fref, 'r') as fh: 
            for cnt, line in enumerate(fh):  # read the input file line by line
                line = line.rstrip().lstrip()
                if refpoint == 'first':
                    if not got_first_point and current_section_type == 'soma':
                        sdg = re_pts.match(line)
                        if sdg is not None:
                            sd = sdg.groups()
                            xyzr = np.array([float(sd[1]), float(sd[2]), float(sd[3]), float(sd[4])])
                            got_first_point = True
                            return xyzr
                else:
                    if current_section_type == 'soma':  # get the LAST soma point
                        sdg = re_pts.match(line)
                        if sdg is not None:
                            sd = sdg.groups()
                            xyzr = np.array([float(sd[1]), float(sd[2]), float(sd[3]), float(sd[4])])
                if current_section_type != 'soma' and xyzr is not None:
                    return xyzr
                sdg = re_append.match(line)
                if sdg is not None:
                    current_section_type = sdg.groups()[0]
                    # print(f"Section type: {current_section_type:s}")
        return None
        
    def print_sections(self, fn=None, fref=None, refpoint='first'):
        if fn is not None:
            self.fn = fn
        if fref is not None:
            refxyzr = self.get_ref_point(fref, refpoint=refpoint)
            print('Base point in reference file (old): ', refxyzr)
        else:
            print('Need reference point')
            exit()
        base_xyz = None
        current_section_type = None
        got_first_point = False
        olines = []
        print('self.fn: ', self.fn)
        with open(self.fn, 'r') as fh: 
            for cnt, line in enumerate(fh):  # read the input file line by line
                line = line.rstrip().lstrip()
                if not got_first_point and current_section_type == 'soma':
                    sdg = re_pts.match(line)
                    if sdg is not None:
                        sd = sdg.groups()
                        xyzr = np.array([float(sd[1]), float(sd[2]), float(sd[3]), float(sd[4])])
                        delta = xyzr - refxyzr
                        got_first_point = True
                        print('base point: ', xyzr)
                        print('delta: ', delta[0:3])
                    olines.append(line)
                elif got_first_point and current_section_type != 'soma':
                    sdg = re_pts.match(line)
                    if sdg is not None:
                        sd = sdg.groups()
                        pt3d = np.array([float(sd[1]), float(sd[2]), float(sd[3]), float(sd[4])])
                        # print(f"{'Old line was':>16s} {current_section_type:>12s}", end='')
                        # print(f" pt3dadd({pt3d[0]:f}, {pt3d[1]:f}, {pt3d[2]:f}, {float(pt3d[3]):f})")
                        pt3d[0:3] += delta[0:3]
                        # print(f"{'New line is':>16s} {' ':>12s} pt3dadd({pt3d[0]:f}, {pt3d[1]:f}, {pt3d[2]:f}, {float(pt3d[3]):f})")
                        olines.append(f"  pt3dadd({pt3d[0]:f}, {pt3d[1]:f}, {pt3d[2]:f}, {float(pt3d[3]):f})")
                    else:
                        olines.append(line)                #

                else:
                    olines.append(line)

                sdg = re_append.match(line)
                if sdg is not None:
                    current_section_type = sdg.groups()[0]
#                    print(f"Section type: {current_section_type:s}"                
                # sdg = re_section.match(line)
                # if sdg is not None:
                #     sd = sdg.groups()
                #     print(f'   sections[{int(sd[1]):3d}]')

        print('Translate Cell, output file: ', 'test.hoc')
        header = f"// Original file: {str(self.fn):s}\n"
        header += '// Coordinates Translated '

        now = datetime.datetime.now()
        header += '\n// on: ' + now.replace(microsecond=0).isoformat(' ') + '\n'
        fno = 'test.hoc'
        with open(fno, 'w') as fh:
            fh.write(header)
            for o in olines:
                fh.write(o+'\n')        
                      
    def translate_cell(self, trans_xyz, fn=None, soma_section=None, reorder=False):
        """
        Translate a cell position by modifying the hoc code
        Writes out a new hoc file with "_translated" in the name
        
        Parameters
        ----------
        trans_xyz: 3 element list or tuple
            Values to translate the position of EVERY point by,
            in order of (x, y, z)
        
        reorder : bool (default: False)
            Cause values on soma section to be reordered in reverse
            order. Used to bring new soma reconstructions in VCNmodel (SBEM data)
            in register with older reconstructions so that connections
            to denrites and axon are correctly positioned.
        
        """
        if fn is not None:
            self.fn = fn
        soma_data = []
        sflag = False
        insoma = False
        soma_top = None
        soma_bottom = None
        fno = Path(self.fn)
        fno_name = fno.stem
        fno_name += '_translated'
        # print('fno: ', fno_name)
        fno = Path(self.fn.parent, fno_name).with_suffix('.hoc')
        print('Translate_cell input file: ', self.fn)
        # print('output fn: ', fno)
        header = f"// Original file: {str(self.fn):s}\n"
        header += '// Coordinates Translated '
        if reorder is False:
            header += 'and in original order'
        else:
            header += 'and soma sections have been reordered'
        now = datetime.datetime.now()
        header += '\n// on: ' + now.replace(microsecond=0).isoformat(' ') + '\n'
        olines = []
        with open(self.fn, 'r') as fh: 
            for cnt, line in enumerate(fh):  # read the input file line by line
                line = line.rstrip().lstrip()
                if line.startswith('soma.append()'):  # beginning of a soma section
                    sflag = True
                # if line in ['hillock.append()', 'axon.append()', 'dendrite.append()',
                #             'unmeyleinatedaxon.append()', 'muyelinatedaxon.append()']:
                #     sflag = False
                if sflag and line.startswith('sections['):  # section definitions will follow
                    insoma = True
                    lcount = 0
                    olines.append(line)
                    somalines = []
                if insoma:
                    sdg = re_pts.match(line)
                    if sdg is not None:
                        sd = sdg.groups()
                        xyzr = [float(sd[1]), float(sd[2]), float(sd[3]), float(sd[4])]
                        soma_data.append(xyzr)
                        if soma_top == None:
                            soma_top = soma_data[-1]
                        xyzr[0:3] -= trans_xyz
                        somalines.append(f"  pt3dadd({xyzr[0]:f}, {xyzr[1]:f}, {xyzr[2]:f}, {float(sd[4]):f})")
                if insoma and line.find('}') != -1:  # end of the section
                    insoma=False
                    if reorder:
                        revsoma = list(reversed(somalines))
                        # print('somalines reordered: ')
                        for o in revsoma:
                            olines.append(o)
                            # print(o)
                    else:
                        # print('somalines normal order: ')
                        for o in somalines:
                            olines.append(o)
                            # print(o)

                if not insoma:
                    olines.append(line)
        # print('\ninput fn: ', self.fn)
        print('Translate Cell, output file: ', str(fno))
        with open(fno, 'w') as fh:
            fh.write(header)
            for o in olines:
                fh.write(o+'\n')
        
    def get_soma(self, fn=None):
        if fn is not None:
            self.fn = fn
        soma_data = []
        sflag = False
        insoma = False
        soma_top = None
        soma_bottom = None
        with open(self.fn, 'r') as fh:
            for cnt, line in enumerate(fh):
                line = line.rstrip().lstrip()
                if line.startswith('soma.append()'):
                    sflag = True
                    continue
                if line in ['hillock.append()', 'axon.append()', 'dendrite.append()', 
                            'unmeyleinatedaxon.append()', 'muyelinatedaxon.append()']:
                            sflag = False
                            continue
                if sflag and line.startswith('sections['):
                    insoma = True
                    lcount = 0
                    
                    continue
                if insoma:
                    sdg = re_pts.match(line)
                    if sdg is not None:
                        sd = sdg.groups()
                        xyzr = [float(sd[1]), float(sd[2]), float(sd[3]), float(sd[4])]
                        soma_data.append(xyzr)
                        if soma_top == None:
                            soma_top = soma_data[-1]
                if insoma and line.find('}') != -1:
                    insoma=False
                    if soma_bottom == None:
                        soma_bottom = soma_data[-1]
        print('\nGet soma from file: ', self.fn)
        # print('soma data: ')
        # for s in soma_data:
        #     print('   ', s)
        # print(soma_top, soma_bottom)
        self.soma_tb = [soma_top, soma_bottom]
        print('soma tb: ', self.soma_tb)
        return(self.soma_tb)
        
        


def main():
    basedir = '/Users/pbmanis/Desktop/Python/VCNModel'
    fn = Path(basedir, 'VCN_Cells', 'VCN_c20', 'Morphology', 'VCN_c20.hoc')
    fn = '/Users/pbmanis/Desktop/Python/VCNModel/VCN_newcoord/VCN_c09.hoc'
    fn = '/Users/pbmanis/Desktop/Python/VCNModel/ASA/CellBodySWCs/VCN_c09_CellBody01.hoc'
    h = ParseHoc(fn)
    h.read_hoc()
    somas = h.get_soma()
    print('somas: ', somas)
if __name__ == '__main__':
    main()
    