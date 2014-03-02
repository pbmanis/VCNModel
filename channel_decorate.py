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
import nrnlibrary.nrnutils as nu
from channel_manager import channelManager

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
        self.biophys(hf )
        hf.update() # make sure we update the information about mechanisms in each section.

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
            print 'Inserting channels as if cell type is %s with modeltype %s ' % (celltype, self.channelInfo.modeltype)
        for s in hf.sections.keys():
            sectype = string.rsplit(s, '[')[0]
            for mech in self.cMan.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    print 'mech %s not found? ' % mech
                    continue
                if verify:
                    print 'section: %s  insert mechanism: %s at ' % (s, mech), self.cMan.channelMap[sectype][mech]
                x = nu.Mechanism(mech)
                x.insert_into(hf.sections[s])
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                gbar = self.gbarAdjust(sectype, mech, s)  # map density by location/distance
                #print 'parmap: ', parMap, mech
                if parMap is not None and mech in parMap.keys():
                    if verify:
                        print 'parMap[mech]', mech, parMap[mech], gbar,
                    gbar = gbar * parMap[mech]
                    if verify:
                        print '  new gbar: ', gbar
                setattr(hf.sections[s](), setup, gbar)
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
        for s in hf.sections.keys():
            sectype = string.rsplit(s, '[')[0]
            print 'Section: %s' % s
            for mech in self.cMan.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    continue
                if verify:
                    print 'Section: %s  find mechanism: %s at ' % (s, mech), self.cMan.channelMap[sectype][mech]
                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                bar = getattr(hf.sections[s](), setup)
                print '\tmech: %8s  gbar: %f' % (mech, bar)
