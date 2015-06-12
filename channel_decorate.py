__author__ = 'pbmanis'

"""
channel_decorate:
Insert biophysical mechanisms into a model.
This function attempts to automatically decorate a hoc-imported model set of sections
with appropriate conductances.

The class takes as input the object hf, which is an instance of the class hoc_reader
It also takes the celltype, a string that directs how conductances should be inserted.

"""
from pylibrary.Params import Params
import pylibrary.PyNrnUtilities as pn
import pprint
import string
import numpy as np
from nrnlibrary.protocols import protocol
from nrnlibrary import cells
from nrnlibrary import synapses
from nrnlibrary.util import get_anspikes
from nrnlibrary.util import sound
import nrnlibrary.util as nu

#import nrnlibrary.nrnutils as nu
from channel_manager import channelManager
excludeMechs = [] # ['ihvcn', 'kht', 'klt', 'nav11']

class ChannelDecorate():
    def __init__(self, hf, celltype=None, modeltype=None, parMap=None):

        self.channelInfo = Params(newCm=1.0,
                              newRa=100.0,  # // changed 10/20/2007 to center in range')
                              newg_leak=0.000004935,
                              eK_def=-85, eNa_def=50,
                              ca_init=70e-6,  # free calcium in molar
                              v_init=-80,  # mV
                              pharmManip={'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False,
                                          'XE': False},
                              celltype=celltype,
                              modeltype=modeltype,
                              distanceMap=hf.distanceMap,
                              parMap=parMap,
        )

        self.cMan = channelManager(celltype+'_'+modeltype)
        self.channelMap = self.cMan.channelMap
        self.distMap = self.cMan.distMap
        self.irange = self.cMan.irange
        #cm.printMap()
        # gmapper allows us to map the names of mechanisms and thier conductance names, which may
        # vary in the hoc files.
        # The versions in the mechanisms directory here have been systematized, but this
        # dictionary may help when adding other conductances.

        self.gmapper = self.cMan.gmapper
        self.biophys(hf)
            

    def biophys(self, hf, verify=False):
        """
        Inputs: run parameter structure, model parameter structure
        verify = True to run through the inserted mechanisms and see that they are really there.
        Outputs: None
        Action: Channel insertion into model
        Side Effects:
            Sets conductances in every different kind of section
            Does not update any class variables (via self).

        original hoc code: Paul B. Manis, Ph.D.
        25 Sept. 2007
        Modified to use gca for HH formulation of calcium current
        14 Oct 2007
        converted for Python, 17 Oct 2012 (PB Manis)
        modified to use new hf hoc_reader class to access section types and mechanisms 10-14 Feb 2014 pbmanis
        """
       # createFlag = False
        celltype = self.channelInfo.celltype
        parMap = self.channelInfo.parMap
        dmap = self.channelInfo.distanceMap
        if self.channelInfo is None:
            raise Exception('biophys - no parameters or info passed!')
        if verify:
            print 'Biophys: Inserting channels as if cell type is %s with modeltype %s ' % (celltype, self.channelInfo.modeltype)
        # print dir(hf)
        # print 'hf section groupkeys: ', hf.sec_groups.keys()
        # print 'hf section groups: ', hf.sec_groups
        
        for s in hf.sec_groups.keys():
            sectype = self.remapSectionType(string.rsplit(s, '[')[0])
            if sectype not in self.cMan.channelMap.keys():
                print 'encountered unknown section group type: %s  Not decorating' % sectype
                continue
            # print 'Biophys: Section type: ', sectype, 'from: ', s
            # print sectype
            # print 'channel mapping keys: ', self.cMan.channelMap.keys()
            for mech in self.cMan.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    print 'mech %s not found? ' % mech
                    continue
                if mech in excludeMechs:
                    continue
                # if mech in ['ihvcn']:
                #     print 'mechanism %s excluded' % mech
                #     continue
                if verify:
                    print 'Biophys: section: %s  insert mechanism: %s at ' % (s, mech), self.cMan.channelMap[sectype][mech]
                x = nu.Mechanism(mech)
                for sec in hf.sec_groups[s]:
                    x.insert_into(hf.get_section(sec))
                    gbar = self.gbarAdjust(sectype, mech, sec)  # map density by location/distance
                #print 'parmap: ', parMap, mech
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                if parMap is not None and mech in parMap.keys():  # note, this allows parmap to have elements BESIDES mechanisms
                    if verify:
                        print 'parMap[mech]', mech, parMap[mech], gbar,
                    gbar = gbar * parMap[mech]
                    if verify:
                        print '  new gbar: ', gbar
                for sec in hf.sec_groups[s]:
                    setattr(hf.get_section(sec), setup, gbar)
        if verify:
            self.channelValidate(hf)


    def gbarAdjust(self, sectype, mech, sec):
        gbar = self.cMan.channelMap[sectype][mech]
        if sectype not in self.cMan.distMap.keys():  # no map for this section type
            return gbar
        elif mech not in self.cMan.distMap[sectype].keys():
            return gbar
        # mecanism exists in the distMap, so we will map gbar to distance from soma
        method = self.cMan.distMap[sectype][mech]['gradient'] # grab the type
        gminf = self.cMan.distMap[sectype][mech]['gminf']
        rate = self.cMan.distMap[sectype][mech]['lambda']
        if method == 'flat':
            return gbar
        if sec in self.channelInfo.distanceMap.keys():
            dist = self.channelInfo.distanceMap[sec]
        else:  # the sec should be in the map, but there could be a coding error that would break that relationship
            raise NameError('gbarAdjust:channel_decorate.py: section %s not in distance map' % sec)
        if method == 'linear':  # rate is "half" point drop
            gbar = gbar - dist*(gbar-gminf)/(2*rate)
            if gbar < 0.:
                gbar = 0. # clip
        elif method in ['exp', 'expdown']:
            gbar = (gbar - gminf) * np.exp(-dist/rate) + gminf
        if gbar < 0.:
            gbar = 0.
        return gbar

    def channelValidate(self, hf, verify=False):
        """
         verify mechanisms insertions -
         go through all the groups, and find inserted conductances and their values
         print the results to the terminal
        """
        secstuff = {}
        for s in hf.sec_groups.keys():
            sectype = self.remapSectionType(string.rsplit(s, '[')[0])
            if sectype not in self.cMan.channelMap.keys():
                print 'encountered unknown section group type: %s  Cannot Validate' % sectype
                continue
#            print 'Validating Section: %s' % s
            for mech in self.cMan.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    continue
                if mech in excludeMechs:
                    continue
                if verify:
                    print '\tSection: %s  find mechanism: %s at ' % (s, mech), self.cMan.channelMap[sectype][mech]
                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                for sec in hf.sec_groups[s]:
                    bar = getattr(hf.get_section(sec), setup)
                    if sec in secstuff.keys():
                        secstuff[sec] += ', g_%s = %g ' % (mech, bar)
                    else:
                        secstuff[sec] = '(%10s) g_%-6s = %g ' % (sectype, mech, bar)
        for i, k in enumerate(secstuff.keys()):
            print '%-20s ' % k, secstuff[k]
                

    def remapSectionType(self, sectype):
        if sectype in ['AXON_0']:
            sectype = 'axon'
        if sectype in ['dendscaled_0', 'dendscaled_1', 'dendscaled_2', 'dendrite']:
            sectype = 'dend'
        return sectype

