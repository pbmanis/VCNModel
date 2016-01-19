__author__ = 'pbmanis'

"""
channel_decorate:
Insert biophysical mechanisms into a model.
This function attempts to automatically decorate a hoc-imported model set of sections
with appropriate conductances.

The class takes as input the object hf, which is an instance of the class hoc_reader
It also takes the celltype, a string that directs how conductances should be inserted.

"""

import string
import numpy as np
from cnmodel.protocols import protocol
from cnmodel import cells
from cnmodel import synapses
from cnmodel.util import get_anspikes
from cnmodel.util import sound
import cnmodel.util as nu

from pylibrary.Params import Params
from cnmodel.util import pynrnutilities
import pprint

"""
chManager is a class to manage channels and their distributions.

For each channel, given a canonical cell type, we define:
basic conductance set
"""
from pylibrary.Params import Params
# import pylibrary.PyNrnUtilities as pn
import pprint as pp
import numpy as np


class ChannelManager():
    def __init__(self, celltype):

        # Scaled params normalized to Rothman and Manis densities:
        self.c_m = 1.0E-6  # in units of F/cm^2
        self.celltype = celltype  # store called cell type.
        if celltype == 'Bushy_RM03':
            totcap = 12.0E-12  # in units of F, from Rothman and Manis, 2003.
            refarea = totcap / self.c_m  # area is in cm^2
            # bushy Rothman-Manis, guinea pig type II
            # model gave cell conductance in nS, but we want S/cm^2 for NEURON
            # so conversion is 1e-9*nS = uS, and refarea is already in cm2
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.0E-9/refarea,
                               kltbar=200.0E-9/refarea,
                               ihbar=20.0E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': self.gBar.kltbar / 4., 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2., 'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar * 0.5, 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar * 0.2, 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4., 'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-1., 1., 11)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nacn': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nacn': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }


        elif celltype in ['TStellate_RM03']:
            totcap = 12.0E-12
            refarea = totcap / self.c_m  # see above for units
            # Type I stellate Rothman and Manis, 2003c
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=0.5E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar / 2.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2, 'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.1, 0.1, 7)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }


        elif celltype in ['TStellate_XM13']:
            totcap = 25.0E-12
            refarea = totcap / self.c_m  # see above for units
            # Type I stellate Rothman and Manis, 2003c
            self.gBar = Params(nabar=800.0E-9/refarea,
                               khtbar=250.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               ihbar=18.0E-9/refarea,
                               leakbar=8.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nav11': 0.0, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nav11': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2, 'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.5, 0.5, 9)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                }



        elif celltype == 'Bushy_XM13':
            # bushy form Xie and Manis, 2013, based on Cao and Oertel mouse conductances
            totcap = 26.0E-12 # uF/cm2 
            refarea = totcap  / self.c_m  # see above for units
            self.gBar = Params(nabar=500.E-9/refarea,
                               khtbar=58.0E-9/refarea,
                               kltbar=80.0E-9/refarea,  # note doubled here... 
                               ihbar=0.25*30.0E-9/refarea,
                               leakbar=0.5*2.0E-9/refarea,
            )
            print 'XM13 gbar:\n', self.gBar.show()
            self.channelMap = {
                'axon': {'nav11': self.gBar.nabar*0, 'klt': self.gBar.kltbar * 0.25, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar * 0.25},
                'hillock': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar*3, 'klt': self.gBar.kltbar*2, 'kht': self.gBar.khtbar*2,
                            'ihvcn': self.gBar.ihbar * 0.5, 'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar * 0.5, 'klt': self.gBar.kltbar *0.5, 'kht': self.gBar.khtbar *0.5,
                         'ihvcn': self.gBar.ihbar *0.5, 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nav11': self.gBar.nabar * 0.25, 'klt': self.gBar.kltbar * 0.25, 'kht': self.gBar.khtbar * 0.25,
                         'ihvcn': self.gBar.ihbar *0.25, 'leak': self.gBar.leakbar * 0.25, },
            }
            self.irange = np.linspace(-2, 2, 7)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.},
                                     'kht': {'gradient': 'llinear', 'gminf': 0., 'lambda': 200.},
                                     'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 200.}}, # gradients are: flat, linear, exponential
                            }

        elif celltype == 'Bushy_XM13PasDend':
            # bushy form Xie and Manis, 2013, based on Cao and Oertel mouse conductances
            # passive dendritestotcap = 26.0E-12 # uF/cm2 
            totcap = 26.0E-12 # uF/cm2 
            refarea = totcap  / self.c_m  # see above for units
            self.gBar = Params(nabar=500.E-9/refarea,
                               khtbar=58.0E-9/refarea,
                               kltbar=80.0E-9/refarea,  # note doubled here... 
                               ihbar=30.0E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            print 'gbar: ', self.gBar
            self.channelMap = {
                'axon': {'nav11': self.gBar.nabar*0, 'klt': self.gBar.kltbar * 0.25, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar * 0.25},
                'hillock': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar*3, 'klt': self.gBar.kltbar*2, 'kht': self.gBar.khtbar*2,
                            'ihvcn': self.gBar.ihbar * 0.5, 'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar * 0.0, 'klt': self.gBar.kltbar*0 , 'kht': self.gBar.khtbar*0,
                         'ihvcn': self.gBar.ihbar*0, 'leak': self.gBar.leakbar*0.5, },
                'apic': {'nav11': self.gBar.nabar * 0.0, 'klt': self.gBar.kltbar * 0, 'kht': self.gBar.khtbar * 0.,
                         'ihvcn': self.gBar.ihbar *0., 'leak': self.gBar.leakbar * 0.25, },
            }
            self.irange = np.linspace(-1, 1, 21)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.},
                                     'kht': {'gradient': 'llinear', 'gminf': 0., 'lambda': 200.},
                                     'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 200.}}, # gradients are: flat, linear, exponential
                            }
        elif celltype == 'Calyx':
            totcap = 12.0E-12
            refarea = totcap / self.c_m  # See above for units # bushy Rothman-Manis, guinea pig type II
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.0E-9/refarea,
                               kltbar=200.0E-9/refarea,
                               ihbar=20.0E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': self.gBar.kltbar / 4., 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2., 'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar * 0.5, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nacn': 0.0, 'klt': self.gBar.kltbar * 0.5, 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': 0.0, 'klt': self.gBar.kltbar * 0.2, 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4., 'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.5, 0.5, 7)
            self.distMap = {}  # uniform, as defined above.


        elif celltype == 'MNTB':
            totcap = 12.0E-12
            refarea = totcap  / self.c_m  # See above for units  # bushy Rothman-Manis, guinea pig type II
            self.gBar = Params(nabar=1000.0E-9/refarea,
                               khtbar=150.00E-9/refarea,
                               kltbar=200.00E-9/refarea,
                               ihbar=20.00E-9/refarea,
                               leakbar=2.00E-9/refarea,
            )
            self.channelMap = {
                  'axon': {'nacn': 0.0, 'klt': self.gBar.kltbar / 4., 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                           'leak': self.gBar.leakbar / 4.},
                  'hillock': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                              'ihvcn': 0., 'leak': self.gBar.leakbar, },
                  'initseg': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                              'ihvcn': self.gBar.ihbar / 2., 'leak': self.gBar.leakbar, },
                  'soma': {'nacn': self.gBar.nabar * 0.5, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                           'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                  'dend': {'nacn': 0.0, 'klt': self.gBar.kltbar * 0.5, 'kht': self.gBar.khtbar * 0.5,
                           'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                  'apic': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar * 0.2, 'kht': self.gBar.khtbar * 0.2,
                           'ihvcn': self.gBar.ihbar / 4., 'leak': self.gBar.leakbar * 0.2, },
              },
            self.irange = np.linspace(-0.3, 0.3, 7)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }

        elif celltype == 'L23pyr':
            self.gBar = Params(nabar=0.30, khtbar=0.002, kltbar=0.0, ihbar=0.0, leakbar=0.0002)
            self.channelMap = {
                'axon': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0., 'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar / 150., 'klt': 0., 'kht': self.gBar.khtbar / 10., 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar / 150., 'klt': 0., 'kht': self.gBar.khtbar / 100., 'ihvcn': self.gBar.ihbar / 3.,
                         'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': self.gBar.nabar / 150, 'klt': 0., 'kht': self.gBar.khtbar / 100., 'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
                }
            self.irange = np.linspace(-0.1, 0.1, 5)
            self.distMap = {}  # uniform

        else:
            raise Exception('Unrecognized cell/parameter set type: %s' % celltype)

        self.gmapper = {'nacn': 'gbar', 'kht': 'gbar', 'klt': 'gbar', 'leak': 'gbar',
                        'ihvcn': 'gbar', 'jsrna': 'gbar', 'nav11': 'gbar'}


    def printMap(self):
        pp.pprint(self.channelMap)



excludeMechs = [] # ['ihvcn', 'kht', 'klt', 'nav11']

class ChannelDecorate():
    def __init__(self, hf, celltype=None, modeltype=None, parMap=None, verify=False):

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

        self.cMan = ChannelManager(celltype+'_'+modeltype)
        self.channelMap = self.cMan.channelMap
        self.distMap = self.cMan.distMap
        self.irange = self.cMan.irange
        #cm.printMap()
        # gmapper allows us to map the names of mechanisms and thier conductance names, which may
        # vary in the hoc files.
        # The versions in the mechanisms directory here have been systematized, but this
        # dictionary may help when adding other conductances.

        self.gmapper = self.cMan.gmapper
        self.hf = self._biophys(hf, verify=verify)
        print 'ChannelDecorate: Model Decorated with channels (if this appears more than once per cell, there is a problem)'


    def _biophys(self, hf, verify=False):
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
            print('Biophys: Inserting channels as if cell type is {:s} with modeltype {:s}'
                 .format(celltype, self.channelInfo.modeltype))
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
                    print('Biophys: section group: {:s}  insert mechanism: {:s} at {:.8f}'
                        .format(s, mech, self.cMan.channelMap[sectype][mech]))
                x = nu.Mechanism(mech)
                for sec in hf.sec_groups[s]:
                    x.insert_into(hf.get_section(sec))
                    if verify:
                        print '   inserting into section', sec
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
        return hf


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
            raise NameError('gbarAdjust in channel_decorate.py: section %s not in distance map' % sec)
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
        print '\nChannel Validation'
        secstuff = {}
        for s in hf.sec_groups.keys():
            sectype = self.remapSectionType(string.rsplit(s, '[')[0])
            if sectype not in self.cMan.channelMap.keys():
                print 'Validation: encountered unknown section group type: %s  Cannot Validate' % sectype
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
#        for i, k in enumerate(secstuff.keys()):
#            print '%-20s ' % k, secstuff[k]
                

    def remapSectionType(self, sectype):
        if sectype in ['AXON_0']:
            sectype = 'axon'
        if sectype in ['dendscaled_0', 'dendscaled_1', 'dendscaled_2', 'dendrite']:
            sectype = 'dend'
        if sectype in ['apical_dendrite']:
            sectype = 'apic'
        return sectype

