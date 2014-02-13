__author__ = 'pbmanis'

"""
channel_decorate:
Insert biophysical mechanisms into a model.
The routine takes as input a dictionary constructed as follows:
{structname: {mech: [name, gbarname]}
"""
from pylibrary.Params import Params
import pylibrary.PyNrnUtilities as pn
import pprint
import string
import nrnlibrary.nrnutils as nu

class ChannelDecorate():
    def __init__(self, cellType, hf, sectionList):

        self.channelInfo = Params(newCm=1.0,
                              newRa=100.0,  # // changed 10/20/2007 to center in range')
                              newg_leak=0.000004935,
                              eK_def=-85, eNa_def=50,
                              ca_init=70e-6,  # free calcium in molar
                              v_init=-80,  # mV
                              pharmManip={'TTX': False, 'ZD': False, 'Cd': False, 'DTX': False, 'TEA': False,
                                          'XE': False},

                              )

        # normalized to Rothman and Manis densities:
        c_m = 1.0
        totcap = 12.0
        refarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
        bg = Params(nabar = pn.nstomho(1000., refarea),
                           khtbar = pn.nstomho(58.0, refarea),
                           kltbar = pn.nstomho(80.0, refarea),
                           ihbar = pn.nstomho(30.0, refarea),
                           leakbar = pn.nstomho(2.0, refarea),
        )
        self.channelMap = {'Bushy': {
            'axon': {'nav11': 0, 'klt': bg.kltbar/2., 'kht': bg.khtbar, 'ihvcn': 0., 'leak': bg.leakbar/4.},
            'hillock': {'nav11': bg.nabar, 'klt': bg.kltbar, 'kht': bg.khtbar, 'ihvcn': 0., 'leak': bg.leakbar,},
            'initseg': {'nav11': bg.nabar, 'klt': bg.kltbar, 'kht': bg.khtbar, 'ihvcn': bg.ihbar/2., 'leak': bg.leakbar,},
            'soma': {'nav11': bg.nabar, 'klt': bg.kltbar, 'kht': bg.khtbar, 'ihvcn': bg.ihbar, 'leak': bg.leakbar,},
            'dend': {'nav11': 0.0, 'klt': bg.kltbar*0.3, 'kht': bg.khtbar*0.1, 'ihvcn': bg.ihbar/3., 'leak': bg.leakbar,},
            'apic': {'nav11': 0.0, 'klt': bg.kltbar*0.1, 'kht': bg.khtbar*0.1, 'ihvcn': bg.ihbar/4., 'leak': bg.leakbar*10.,},

            },
        }
        # self.channelMap = {'Bushy': {
        #    'axon': {'nav11': bg.nabar,}, #'leak': bg.leakbar/4.}, #'nav11': 0., 'klt': 01., 'kht': bg.khtbar, 'ihvcn': 0., },
        #    'hillock': {'nav11': bg.nabar,}, #'leak': bg.leakbar,}, #'nav11': bg.nabar, 'klt': bg.kltbar, 'kht': bg.khtbar, 'ihvcn': 0., 'leak': bg.leakbar,},
        #    'initseg': {'nav11': bg.nabar,}, #'leak': bg.leakbar,}, #'nav11': bg.nabar, 'klt': bg.kltbar, 'kht': bg.khtbar, 'ihvcn': bg.ihbar/2., 'leak': bg.leakbar,},
        #    'soma': {'nav11': bg.nabar,}, #{'leak': bg.leakbar*10,}, # 'ih': bg.ihbar}, #'klt': bg.kltbar, }, #'kht': bg.khtbar,  #  'inav11': bg.nabar, }},
        #    'dend': {'nav11': bg.nabar,}, #'leak': bg.leakbar,}, #'nav11': 0.0, 'klt': bg.kltbar*0.3, 'kht': bg.khtbar*0.1, 'ihvcn': bg.ihbar/3., 'leak': bg.leakbar,},
        #    'apic': {'nav11': bg.nabar,}, #'leak': bg.leakbar,}, #'nav11': 0.0, 'klt': bg.kltbar*0.1, 'kht': bg.khtbar*0.1, 'ihvcn': bg.ihbar/4., 'leak': bg.leakbar*10.,},
        #
        #     },
        # }

        gmapper = {'nav11': 'gbar', 'kht': 'gbar', 'klt': 'gbar', 'leak': 'gbar', 'ihvcn': 'gbar'}
        self.biophys(hf, cellType, gmapper)
        hf.update() # make sure we update the information about mechanisms in each section.

    def biophys(self, hf, cellType, gmapper, createFlag=True):
        """
        Inputs: run parameter structure, model parameter structure
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
        """
       # createFlag = False
        if self.channelInfo is None:
            raise Exception('biophys - no parameters or info passed!')
        for s in hf.sections.keys():
            sectype = string.rsplit(s, '[')[0]
            for mech in self.channelMap[cellType][sectype].keys():
                if mech not in gmapper.keys():
                    continue
                x = nu.Mechanism(mech) # , {gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                x.set_parameters({gmapper[mech]: self.channelMap[cellType][sectype][mech]})
                x.insert_into(hf.sections[s])


        # if createFlag:
        #     for part in hf.sec_groups.keys():
        #         if part not in self.channelMap[cellType].keys():
        #             continue
        #         for sectionID in hf.sec_groups[part]:
        #             sec = eval('hf.h.%s' % (sectionID))
        #             for mech in self.channelMap[cellType][part].keys():
        #                 if mech not in gmapper.keys():
        #                     continue
        #                 try:
        #                     sec.insert(mech)
        #                 except:
        #                     print 'missing mech: ', mech
        #                 setup = ('%s_%s' % (gmapper[mech], mech))
        #                 setattr(sec(),setup, self.channelMap[cellType][part][mech])


