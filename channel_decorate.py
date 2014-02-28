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
    def __init__(self, hf, celltype=None, modeltype=None):

        self.channelInfo = Params(newCm=1.0,
                              newRa=100.0,  # // changed 10/20/2007 to center in range')
                              newg_leak=0.000004935,
                              eK_def=-85, eNa_def=50,
                              ca_init=70e-6,  # free calcium in molar
                              v_init=-80,  # mV
                              pharmManip={'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False,
                                          'XE': False},
                              celltype=celltype,
                              modeltype = modeltype,
        )

        cm = channelManager(celltype+'_'+modeltype)
        self.channelMap = cm.channelMap
        self.irange = cm.irange
        #cm.printMap()
        # gmapper allows us to map the names of mechanisms and thier conductance names, which may
        # vary in the hoc files.
        # The versions in the mechanisms directory here have been systematized, but this
        # dictionary may help when adding other conductances.

        self.gmapper = cm.gmapper
        self.biophys(hf)
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
        if self.channelInfo is None:
            raise Exception('biophys - no parameters or info passed!')
        if verify:
            print 'Inserting channels as if cell type is %s with modeltype %s ' % (celltype, self.channelInfo.modeltype)
        for s in hf.sections.keys():
            sectype = string.rsplit(s, '[')[0]
            for mech in self.channelMap[sectype].keys():
                if mech not in self.gmapper.keys():
                    print 'mech %s not found? ' % mech
                    continue
                if verify:
                    print 'section: %s  insert mechanism: %s at ' % (s, mech), self.channelMap[sectype][mech]

                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                x.insert_into(hf.sections[s])
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                setattr(hf.sections[s](),setup, self.channelMap[sectype][mech])

                #x.set_parameters({gmapper[mech]: self.channelMap[celltype][sectype][mech]})

        # # old method : the new method is cleaner.
        # for part in hf.sec_groups.keys():
        #     if part not in self.channelMap[celltype].keys():
        #         continue
        #     for sectionID in hf.sec_groups[part]:
        #         sec = eval('hf.h.%s' % (sectionID))
        #         for mech in self.channelMap[celltype][part].keys():
        #             if mech not in gmapper.keys():
        #                 continue
        #             try:
        #                 sec.insert(mech)
        #                 print 'inserted %s into sec= ' % (mech), part, sectionID
        #             except:
        #                 print 'missing mech: ', mech
        #             setup = ('%s_%s' % (gmapper[mech], mech))
        #             setattr(sec(),setup, self.channelMap[celltype][part][mech])
        if verify:
            self.channelValidate(hf)

    def channelValidate(self, hf):     # verify insertions - go through all the groups, and find inserted conductances and thier values
        for part in hf.sec_groups.keys():
            if part not in self.channelMap[self.channelInfo.celltype].keys():
                continue
            for sectionID in hf.sec_groups[part]:
                sec = eval('hf.h.%s' % (sectionID))
                for mech in self.channelMap[self.channelInfo.celltype][part].keys():
                    if mech not in self.gmapper.keys():
                        continue
                    setup = ('%s_%s' % (self.gmapper[mech], mech))
                    g = getattr(sec(),setup, self.channelMap[self.channelInfo.celltype][part][mech])
                    print '%s %s %s = %g' % (part, sectionID, mech, g)
