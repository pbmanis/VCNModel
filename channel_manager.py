__author__ = 'pbmanis'
"""
chManager is a class to manage channels and their distributions.

For each channel, given a canonical cell type, we define:
basic conductance set
"""
from pylibrary.Params import Params
import pylibrary.PyNrnUtilities as pn
import pprint as pp
import numpy as np


class channelManager():
    def __init__(self, celltype):

        # Scaled params normalized to Rothman and Manis densities:
        self.c_m = 1.0E-6  # in units of F/cm^2

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
            self.irange = np.linspace(-2., 2., 7)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nacn': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'nacn': {'gradient': 'exp', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }


        elif celltype in ['Stellate_RM03']:
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


        elif celltype in ['Stellate_XM13']:
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
                               ihbar=30.0E-9/refarea,
                               leakbar=2.0E-9/refarea,
            )
            print 'gbar: ', self.gBar
            self.channelMap = {
                'axon': {'nav11': self.gBar.nabar*0.5, 'klt': self.gBar.kltbar * 0.25, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                         'leak': self.gBar.leakbar * 0.25},
                'hillock': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar * 0.5, 'leak': self.gBar.leakbar, },
                'soma': {'nav11': self.gBar.nabar, 'klt': self.gBar.kltbar, 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar, },
                'dend': {'nav11': self.gBar.nabar * 0.5, 'klt': self.gBar.kltbar , 'kht': self.gBar.khtbar,
                         'ihvcn': self.gBar.ihbar, 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nav11': self.gBar.nabar * 0.25, 'klt': self.gBar.kltbar * 0.5, 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar *0.25, 'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-1, 1, 7)
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
