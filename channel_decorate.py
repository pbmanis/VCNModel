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
import nrnlibrary.nrnutils as nu

class ChannelDecorate():
    def __init__(self, hf, celltype=None):

        self.channelInfo = Params(newCm=1.0,
                              newRa=100.0,  # // changed 10/20/2007 to center in range')
                              newg_leak=0.000004935,
                              eK_def=-85, eNa_def=50,
                              ca_init=70e-6,  # free calcium in molar
                              v_init=-80,  # mV
                              pharmManip={'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False,
                                          'XE': False},

                              celltype=celltype,
        )

        # Scaled parames normalized to Rothman and Manis densities:
        c_m = 1.0
        totcap = 12.0
        refarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
        # bushy Rothman-Manis, guinea pig type II
        b_rm = Params(nabar = pn.nstomho(1000., refarea),
                           khtbar = pn.nstomho(150.0, refarea),
                           kltbar = pn.nstomho(200.0, refarea),
                           ihbar = pn.nstomho(20.0, refarea),
                           leakbar = pn.nstomho(2.0, refarea),
        )
        # Type I stellate Rothman and Manis, 2003c
        st_rm = Params(nabar = pn.nstomho(1000., refarea),
                           khtbar = pn.nstomho(150.0, refarea),
                           kltbar = pn.nstomho(0.0, refarea),
                           ihbar = pn.nstomho(0.5, refarea),
                           leakbar = pn.nstomho(2.0, refarea),
        )

        # bushy form Xie and Manis, 2013, based on Cao and Oertel mouse conductances
        totcap = 26.0
        refarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
        b_xm = Params(nabar = pn.nstomho(1000., refarea),
                           khtbar = pn.nstomho(58.0, refarea),
                           kltbar = pn.nstomho(80.0, refarea),
                           ihbar = pn.nstomho(30.0, refarea),
                           leakbar = pn.nstomho(8.0, refarea),
        )
        l23_ms = Params(nabar=0.30, khtbar=0.002, kltbar=0.0, ihbar=0.0, leakbar=0.0002)
        self.channelMap = {'Bushy_RM03': {
            'axon': {'nacn': 0.0, 'klt': b_rm.kltbar/4., 'kht': b_rm.khtbar, 'ihvcn': 0., 'leak': b_rm.leakbar/4.},
            'hillock': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': 0., 'leak': b_rm.leakbar,},
            'initseg': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': b_rm.ihbar/2., 'leak': b_rm.leakbar,},
            'soma': {'nacn': b_rm.nabar*0.5, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': b_rm.ihbar, 'leak': b_rm.leakbar,},
            'dend': {'nacn': 0.0, 'klt': b_rm.kltbar*0.5, 'kht': b_rm.khtbar*0.5, 'ihvcn': b_rm.ihbar/3., 'leak': b_rm.leakbar*0.5,},
            'apic': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar*0.2, 'kht': b_rm.khtbar*0.2, 'ihvcn': b_rm.ihbar/4., 'leak': b_rm.leakbar*0.2,},
            },
            'Bushy_XM13': {
            'axon': {'nacn': 0.0, 'klt': b_xm.kltbar/4., 'kht': b_xm.khtbar, 'ihvcn': 0., 'leak': b_xm.leakbar/4.},
            'hillock': {'nacn': b_xm.nabar, 'klt': b_xm.kltbar, 'kht': b_xm.khtbar, 'ihvcn': 0., 'leak': b_xm.leakbar,},
            'initseg': {'nacn': b_xm.nabar, 'klt': b_xm.kltbar, 'kht': b_xm.khtbar, 'ihvcn': b_xm.ihbar/2., 'leak': b_xm.leakbar,},
            'soma': {'nacn': b_xm.nabar*0.5, 'klt': b_xm.kltbar, 'kht': b_xm.khtbar, 'ihvcn': b_xm.ihbar, 'leak': b_xm.leakbar,},
            'dend': {'nacn': 0.0, 'klt': b_xm.kltbar*0.5, 'kht': b_xm.khtbar*0.5, 'ihvcn': b_xm.ihbar/3., 'leak': b_xm.leakbar*0.5,},
            'apic': {'nacn': b_xm.nabar, 'klt': b_xm.kltbar*0.2, 'kht': b_xm.khtbar*0.2, 'ihvcn': b_xm.ihbar/4., 'leak': b_xm.leakbar*0.2,},
            },
            'Calyx': {
            'axon': {'nacn': 0.0, 'klt': b_rm.kltbar/4., 'kht': b_rm.khtbar, 'ihvcn': 0., 'leak': b_rm.leakbar/4.},
            'hillock': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': 0., 'leak': b_rm.leakbar,},
            'initseg': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': b_rm.ihbar/2., 'leak': b_rm.leakbar,},
            'soma': {'nacn': b_rm.nabar*0.5, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': b_rm.ihbar, 'leak': b_rm.leakbar,},
            'dend': {'nacn': 0.0, 'klt': b_rm.kltbar*0.5, 'kht': b_rm.khtbar*0.5, 'ihvcn': b_rm.ihbar/3., 'leak': b_rm.leakbar*0.5,},
            'apic': {'nacn': 0.0, 'klt': b_rm.kltbar*0.2, 'kht': b_rm.khtbar*0.2, 'ihvcn': b_rm.ihbar/4., 'leak': b_rm.leakbar*0.2,},
            },
            'Stellate': {
            'axon': {'nacn': 0.0, 'klt': 0., 'kht': st_rm.khtbar, 'ihvcn': 0., 'leak': st_rm.leakbar/4.},
            'hillock': {'nacn': st_rm.nabar, 'klt': 0., 'kht': st_rm.khtbar, 'ihvcn': 0., 'leak': st_rm.leakbar,},
            'initseg': {'nacn': st_rm.nabar, 'klt': 0., 'kht': st_rm.khtbar, 'ihvcn': st_rm.ihbar/2., 'leak': st_rm.leakbar,},
            'soma': {'nacn': st_rm.nabar*0.5, 'klt': 0., 'kht': st_rm.khtbar, 'ihvcn': st_rm.ihbar, 'leak': st_rm.leakbar,},
            'dend': {'nacn': 0.0, 'klt': st_rm.kltbar*0., 'kht': st_rm.khtbar*0.5, 'ihvcn': st_rm.ihbar/3., 'leak': st_rm.leakbar*0.5,},
            'apic': {'nacn': 0.0, 'klt': st_rm.kltbar*0., 'kht': st_rm.khtbar*0.2, 'ihvcn': st_rm.ihbar/4., 'leak': st_rm.leakbar*0.2,},
            },
           'MNTB': {
            'axon': {'nacn': 0.0, 'klt': b_rm.kltbar/4., 'kht': b_rm.khtbar, 'ihvcn': 0., 'leak': b_rm.leakbar/4.},
            'hillock': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': 0., 'leak': b_rm.leakbar,},
            'initseg': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': b_rm.ihbar/2., 'leak': b_rm.leakbar,},
            'soma': {'nacn': b_rm.nabar*0.5, 'klt': b_rm.kltbar, 'kht': b_rm.khtbar, 'ihvcn': b_rm.ihbar, 'leak': b_rm.leakbar,},
            'dend': {'nacn': 0.0, 'klt': b_rm.kltbar*0.5, 'kht': b_rm.khtbar*0.5, 'ihvcn': b_rm.ihbar/3., 'leak': b_rm.leakbar*0.5,},
            'apic': {'nacn': b_rm.nabar, 'klt': b_rm.kltbar*0.2, 'kht': b_rm.khtbar*0.2, 'ihvcn': b_rm.ihbar/4., 'leak': b_rm.leakbar*0.2,},
            },
           'L23pyr': {
            'axon': {'nacn': l23_ms.nabar, 'klt': 0., 'kht': l23_ms.khtbar, 'ihvcn': 0., 'leak': l23_ms.leakbar/4.},
            'hillock': {'nacn': l23_ms.nabar, 'klt': 0., 'kht': l23_ms.khtbar, 'ihvcn': 0., 'leak': l23_ms.leakbar,},
            'initseg': {'nacn': l23_ms.nabar, 'klt': 0., 'kht': l23_ms.khtbar, 'ihvcn': l23_ms.ihbar/2., 'leak': l23_ms.leakbar,},
            'soma': {'nacn': l23_ms.nabar/150., 'klt': 0., 'kht': l23_ms.khtbar/10., 'ihvcn': l23_ms.ihbar, 'leak': l23_ms.leakbar,},
            'dend': {'nacn': l23_ms.nabar/150., 'klt': 0., 'kht': l23_ms.khtbar/100., 'ihvcn': l23_ms.ihbar/3., 'leak': l23_ms.leakbar*0.5,},
            'apic': {'nacn': l23_ms.nabar/150, 'klt': 0., 'kht': l23_ms.khtbar/100., 'ihvcn': l23_ms.ihbar/4., 'leak': l23_ms.leakbar*0.2,},
            },
           }

        # gmapper allows us to map the names of mechanisms and thier conductance names, which may
        # vary in the hoc files.
        # The versions in the mechanisms directory here have been systematized, but this
        # dictionary may help when adding other conductances.
        self.gmapper = {'nacn': 'gbar', 'kht': 'gbar', 'klt': 'gbar', 'leak': 'gbar', 'ihvcn': 'gbar'}
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
            print 'Inserting channels as if cell type is %s' % (celltype)
        for s in hf.sections.keys():
            sectype = string.rsplit(s, '[')[0]
            for mech in self.channelMap[celltype][sectype].keys():
                if mech not in self.gmapper.keys():
                    continue
                if verify:
                    print 'section: %s  insert mechanism: %s at ' % (s, mech), self.channelMap[celltype][sectype][mech]

                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                x.insert_into(hf.sections[s])
                setup = ('%s_%s' % (self.gmapper[mech], mech))
                setattr(hf.sections[s](),setup, self.channelMap[celltype][sectype][mech])

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
