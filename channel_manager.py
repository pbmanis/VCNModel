__author__ = 'pbmanis'
"""
chManager is a class to manage channels and their distributions.

For each channel, given a canonical cell type, we define:
basic conductance set
"""
from pylibrary.Params import Params
import pylibrary.PyNrnUtilities as pn


class channelManager():
    def __init__(self, celltype):

        # Scaled params normalized to Rothman and Manis densities:
        self.c_m = 1.0

        if celltype == 'Bushy_RM03':
            totcap = 12.0
            refarea = totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
            # bushy Rothman-Manis, guinea pig type II
            self.gBar = Params(nabar=pn.nstomho(['na'], refarea),
                               khtbar=pn.nstomho(150.0, refarea),
                               kltbar=pn.nstomho(200.0, refarea),
                               ihbar=pn.nstomho(20.0, refarea),
                               leakbar=pn.nstomho(2.0, refarea),
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
                'apic': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar * 0.2, 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4., 'leak': self.gBar.leakbar * 0.2, },
            }

        elif celltype == 'Stellate':
            totcap = 12.0
            refarea = totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
            # Type I stellate Rothman and Manis, 2003c
            self.gBar = Params(nabar=pn.nstomho(1000., refarea),
                               khtbar=pn.nstomho(150.0, refarea),
                               kltbar=pn.nstomho(0.0, refarea),
                               ihbar=pn.nstomho(0.5, refarea),
                               leakbar=pn.nstomho(2.0, refarea),
            )
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': 0.,
                            'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar / 2.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2, 'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }

        elif celltype == 'Bushy_XM13':
            # bushy form Xie and Manis, 2013, based on Cao and Oertel mouse conductances
            totcap = 26.0
            refarea = totcap * 1E-6 / c_m  # pf -> uF, cm = 1uf/cm^2 nominal
            self.gBar = Params(nabar=pn.nstomho(1000., refarea),
                               khtbar=pn.nstomho(58.0, refarea),
                               kltbar=pn.nstomho(80.0, refarea),
                               ihbar=pn.nstomho(30.0, refarea),
                               leakbar=pn.nstomho(8.0, refarea),
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
                'apic': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar * 0.2, 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4., 'leak': self.gBar.leakbar * 0.2, },
            }
        elif celltype == 'Calyx':
            totcap = 12.0
            refarea = totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal  # bushy Rothman-Manis, guinea pig type II
            self.gBar = Params(nabar=pn.nstomho(['na'], refarea),
                               khtbar=pn.nstomho(150.0, refarea),
                               kltbar=pn.nstomho(200.0, refarea),
                               ihbar=pn.nstomho(20.0, refarea),
                               leakbar=pn.nstomho(2.0, refarea),
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

        elif celltype == 'MNTB':
            totcap = 12.0
            refarea = totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal  # bushy Rothman-Manis, guinea pig type II
            self.gBar = Params(nabar=pn.nstomho(['na'], refarea),
                               khtbar=pn.nstomho(150.0, refarea),
                               kltbar=pn.nstomho(200.0, refarea),
                               ihbar=pn.nstomho(20.0, refarea),
                               leakbar=pn.nstomho(2.0, refarea),
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

        else:
            raise Exception('Unrecognized cell/parameter set type: %s' % celltype)

