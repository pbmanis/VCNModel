# -*- coding: utf8 -*-
import numpy as np
import os.path
from pathlib import Path

"""
SWC File format from CNIC:

n T x y z R P
n is an integer label that identifies the current point and increments by one from one line to the next.

T is an integer representing the type of neuronal segment, such as soma, axon, apical dendrite, etc. The standard accepted integer values are given below.

0 = undefined
1 = soma
2 = axon
3 = dendrite
4 = apical dendrite
5 = fork point
6 = end point
7 = custom
x, y, z gives the cartesian coordinates of each node.

R is the radius at that node.

P indicates the parent (the integer label) of the current point or -1 to indicate an origin (soma).
"""


class SWC(object):
    """Encapsulates a morphology tree as defined by the SWC standard.
    
    Parameters
    ----------
    filename : str or None
        The name of an swc file to load
    types : dict or None
        A dictionary mapping {type_id: type_name} that describes the type IDs
        in the swc data (second column).
    data : ndarray or None
        Optionally, a data array may be provided instead of an swc file. This
        is used internally.
    """
    def __init__(self, filename=None, types=None, data=None):
        self._dtype = [
            ('id', int), 
            ('type', int), 
            ('x', float), 
            ('y', float), 
            ('z', float), 
            ('r', float), 
            ('parent', int)
        ]
        
        self._id_lookup = None
        self._sections = None
        self._children = None

        self.sectypes = {
            0: 'undefined',
            1: 'soma',
            2: 'axon',
            3: 'basal_dendrite',
            4: 'apical_dendrite',
        }
        if types is not None:
            self.sectypes.update(types)
        
        if data is not None:
            self.data = data
        elif filename is not None:
            self.load(filename)
        else:
            raise TypeError("Must initialize with filename or data array.")
        
        self.sort()
        
    def load(self, filename):
        self.filename = filename
        print('loading: ', filename)
        self.data = np.loadtxt(filename, dtype=self._dtype)

        
    def copy(self):
        return SWC(data=self.data.copy(), types=self.sectypes)

    @property
    def lookup(self):
        """Return a dict that maps *id* to *index* in the data array.
        """
        if self._id_lookup is None:
            self._id_lookup = dict([(rec['id'], i) for i, rec in enumerate(self.data)])
            #self._id_lookup = {}
            #for i, rec in enumerate(self.data):
                #self._id_lookup[rec['id']] = i
        return self._id_lookup

    def children(self, id):
        """Return a list of all children of the node *id*.
        """
        if self._children is None:
            self._children = {}
            for rec in self.data:
                ch = self._children.setdefault(rec['parent'], [])
                ch.append(rec['id'])
        return self._children.get(id, [])

    def __getitem__(self, id):
        """Return record for node *id*.
        """
        return self.data[self.lookup[id]]

    def reparent(self, id):
        """Rearrange tree to make *id* the new root parent.
        """
        d = self.data
        
        # bail out if this is already the root
        if self[id]['parent'] == -1:
            return
        
        parent = -1
        while id != -1:
            oldparent = self[id]['parent']
            self[id]['parent'] = parent
            parent = id
            id = oldparent
            
        self._children = None
        self.sort()
        
    @property
    def sections(self):
        """Return lists of IDs grouped by topological section.
        
        The first item in each list connects to the last item in a previous
        list.
        """
        if self._sections is None:
            sections = []
            sec = []
            
            # find all nodes with nore than 1 child
            branchpts = set()
            endpoints = set(self.data['id'])
            endpoints.add(-1)
            seen = set()
            for r in self.data:
                p = r['parent']
                if p in seen:
                    branchpts.add(p)
                else:
                    seen.add(p)
                    endpoints.remove(p)
            
            # build lists of unbranched node chains
            for r in self.data:
                sec.append(r['id'])
                if r['id'] in branchpts or r['id'] in endpoints:
                    sections.append(sec)
                    sec = []
            
            self._sections = sections
            
        return self._sections
        
    def connect(self, parent_id, swc):
        """Combine this tree with another by attaching the root of *swc* as a 
        child of *parent_id*.
        """
        data = swc.data.copy()
        shift = self.data['id'].max() + 1 - data['id'].min()
        data['id'] += shift
        rootmask = data['parent'] == -1
        data['parent'] += shift
        data['parent'][rootmask] = parent_id
        
        self.data = np.concatenate([self.data, data])
        self._children = None
        self.sort()
        
    def set_type(self, typ):
        self.data['type'] = typ
        
    def write_hoc(self, filename, types=None):
        """Write data to a HOC file.
        
        Each node type is written to a separate section list.
        """
        hoc = []
        
        sectypes = self.sectypes.copy()
        for t in np.unique(self.data['type']):
            if t not in sectypes:
                sectypes[t] = 'type_%d' % t

        # create section lists
        for t in sectypes.values():
            hoc.extend(['objref %s' % t,
                        '%s = new SectionList()' % t])
        hoc.append('')
            
        # create sections
        sects = self.sections
        hoc.append('create sections[%d]' % len(sects))
        sec_ids = {}
        for i, sec in enumerate(sects):
            # remember hoc index for this section
            endpt = self[sec[-1]]['id']
            sec_id = len(sec_ids)
            sec_ids[endpt] = sec_id
            
            # add section to list
            hoc.append('access sections[%d]' % sec_id)
            typ = self[sec[0]]['type']
            hoc.append('%s.append()' % sectypes[typ])
            
            # connect section to parent
            p = self[sec[0]]['parent']
            if p != -1:
                hoc.append('connect sections[%d](0), sections[%d](1)' % (sec_id, sec_ids[p]))

            # set up geometry for this section
            hoc.append('sections[%d] {' % sec_id)
            for seg in sec:
                rec = self[seg]
                hoc.append('  pt3dadd(%f, %f, %f, %f)' % (rec['x'], rec['y'], rec['z'], rec['r']*2))
            hoc.append('}')
            
            hoc.append('')
        
        open(filename, 'w').write('\n'.join(hoc))

    @property
    def root(self):
        """ID of the root node of the tree.
        """
        print(self.data['parent'])
        ind = np.argwhere(self.data['parent'] == -1)[0, 0]
        return self.data[ind]['id']

    def sort(self):
        """Sort the tree in topological order.
        """
        order = self.branch(self.root)
        lt = self.lookup
        indexes = np.array([lt[i] for i in order], dtype=int)
        self.data = self.data[indexes]
        
        self._id_lookup = None
        self._sections = None
        
    def path(self, node):
        path = [node]
        while True:
            node = self[node]['parent']
            if node < 0:
                return path
            path.append(node)

    def scale(self, x, y, z, r):
        self.data['x'] *= x
        self.data['y'] *= y
        self.data['z'] *= z
        self.data['r'] *= r
        
    def translate(self, x, y, z):
        self.data['x'] += x
        self.data['y'] += y
        self.data['z'] += z
        
    def branch(self, id):
        """Return a list of IDs in the branch beginning at *id*.
        """
        branch = [id]
        for ch in self.children(id):
            branch.extend(self.branch(ch))
        return branch
    
    def topology(self):
        """Print the tree topology.
        """
        path = []
        indent = ''
        secparents = [self[s[0]]['parent'] for s in self.sections]
        
        for i, sec in enumerate(self.sections):
            p = secparents[i]
            if p != -1:
                ind = path.index(p)
                path = path[:ind+1]
                indent = indent[:(ind+1) * 3]
            path.append(self[sec[-1]]['id'])

            # look ahead to see whether subsequent sections are children
            if p in secparents[i+1:]:
                this_indent = indent[:-2] + u"├─ "
                indent =      indent[:-2] + u"│  │  "
            else:
                this_indent = indent[:-2] + u"└─ "
                indent =      indent[:-2] + u"   │  "
                
                
            typ = self.sectypes[self[sec[0]]['type']]
            if len(sec) > 10:
                secstr = "%s,...%s" % (str(tuple(sec[:3]))[:-1], str(tuple(sec[-3:]))[1:])
            else:
                secstr = str(tuple(sec))
            print("%ssections[%d] type=%s parent=%d %s" % (this_indent, i, typ, p, secstr))


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 1:
        exit()
    if sys.argv[1] == 'c08':
        base = 'VCN_Cells/VCN_c08/Morphology'
        somafile = os.path.join(base, 'VCN_c08_cellbody.swc')
        soma = SWC(somafile, types={1:'soma', 2:'axon', 3:'dendrite'})
        soma.set_type(1)
        soma.data['r'] *= 0.5  # this data was recorded as diameter
        soma.data['r'] *= 0.54545 # per michael, but not clear if this is really true or not?
     
        axon = SWC(os.path.join(base, 'VCN_c08_axon_hillock.swc'))
        axon.set_type(2)
        dend = SWC(os.path.join(base, 'VCN_c18_dendnonscaled.swc'))
        dend.set_type(3)
        dend.reparent(755)
    
        cell = soma.copy()
        cell.connect(57, axon)
        cell.connect(39, dend)
        # correct for voxel size
        cell.scale(0.11, 0.11, 0.06, 0.11)
        # correct for shrinkage
        s = 1.0 / 0.86 # 0.75
        cell.scale(s, s, s, s)
        cell.translate(-70, -90, -60)
        cell.write_hoc(os.path.join(base, 'VCN_c08.hoc'))
        soma.topology()

        #soma.write_hoc('MorphologyFiles/VCN_c18_reparented755.hoc')
    elif sys.argv[1] == 'P30':
        print ('P30')
        base = 'Calyx_Terminals/VCN_P30/Morphology'
        # somafile = os.path.join(base, 'P30_cellbody_scaled.swc')
        # soma = SWC(somafile, types={1:'soma', 2:'axon', 3:'dendrite'})
        # soma.set_type(1)
        # soma.data['r'] *= 0.1  # this data was recorded as diameter
     
        axon = SWC(os.path.join(base, 'input1_axon_scaled.swc'))
        axon.set_type(2)
        axon.scale(0.11, 0.11, 0.11, 0.11/2.)
        axon.translate(-20, -180, -20)
        collateral = SWC(os.path.join(base, 'input1_collat_scaled.swc'))
        collateral.set_type(2)
        collateral.scale(0.11, 0.11, 0.11, 0.11/2.)
        collateral.translate(-20, -180, -20)
        calyx = SWC(os.path.join(base, 'input_scaled.swc'))
        calyx.scale(0.11, 0.11, 0.11, 0.11/2.)
        calyx.translate(-20, -180, -20)
        calyx.set_type(3)
        #dend.reparent(755)
        calyxbody = calyx.copy()
        calyxbody.connect(1, axon)
        #axon.connect(31, collateral)
#       calyx.connect(78, dend)
        # correct for shrinkage
        #s = 1.0 / 0.86 # 0.75
        #cell.scale(s, s, s, s)
        calyx.write_hoc(os.path.join('MorphologyFiles', 'P30_calyx.hoc'))
        
        calyx.topology()
    else:
        fn = Path(sys.argv[1])
        cell = SWC(str(fn))
        cell.write_hoc(Path(fn.parent, fn.name+'.hoc'))
