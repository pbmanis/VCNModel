"""  Test shuffled autorrelation. (for pytest)
Shuffled autocorrelatioon routine testing We confirm that the status of the
results for all 6 patterns is the same as defined in the audit. The test data
are computed in the sac.py module

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021- Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
import numpy as np
import vcnmodel.analyzers.sac as sac
from vcnmodel.util.user_tester import UserTester

SAC = sac.SAC()


def test_sac():
    SACTester(key="ShuffledAutocorrelation")


class SACTester(UserTester):
    """
    Test the Shuffled Autocorrelation routine. We only test the cython version
    here, but have verifie that the reference versions in Python and Numba yield
    exactly the same results (see sac.py) Note that in sac.makeTests, the random
    number generator is seeded with the same value for every run, which should
    help with consistency in the unit tests.
    """

    def __init__(self, key):
        UserTester.__init__(self, key)

    def run_test(self, plotflag=False):
        """
        Run the tests - doing ALL of the panels shown in
        the Louage et al. 2004 paper, fFgure 2.
        see sac.py for plots of these tests.
        """
        self.plotflag = plotflag
        tests = ["A", "B", "C", "D", "E", "F", "G"]
        ntests = len(tests)
        maxn = 25000000
        yr = [[]] * ntests  # store the SAC results for each test
        xr = [[]] * ntests
        sacr = [[]] * ntests
        for i, t in enumerate(tests):
            X = SAC.makeTests(t)
            yr[i], xr[i] = SAC.SAC_with_histo(X, SAC.SPars, engine="cython", dither=0.0)
            save_indices = np.where((xr[i] >= -0.005) & (xr[i] <= 0.005))
            yr[i] = yr[i][save_indices]
            xr[i] = xr[i][save_indices]
            sacr[i] = SAC.SACResult

        max_sac = np.zeros(ntests)
        min_sac = np.zeros(ntests)
        mean_sac = np.zeros(ntests)
        zero_sac = np.zeros(ntests)
        for i in range(len(xr)):
            max_sac[i] = np.max(yr[i])
            min_sac[i] = np.min(yr[i])
            mean_sac[i] = np.mean(yr[i])
            xzero = int(np.where(xr[i] == 0.0)[0])
            zero_sac[i] = yr[i][xzero]  # save the value at 0 time
        self.yr = yr
        self.xr = xr
        self.tests = tests
        if self.audit and plotflag:
            self.show_result()

        info = dict(
            max_sac=max_sac,
            min_sac=min_sac,
            mean_sac=mean_sac,
            zero_sac=zero_sac,
        )
        return info

    def show_result(self):
        """
        Simple display of the sac results
        Use sac.py for finer-grained plots.
        """
        import matplotlib.pyplot as mpl

        fig, ax = mpl.subplots(len(self.tests), 1)
        ax = ax.ravel()
        for i in range(len(self.tests)):
            ax[i].plot(self.xr[i], self.yr[i], "k-")
            ax[i].set_xlim(-0.005, 0.005)
        mpl.show()
        return

    def assert_test_info(self, *args, **kwds):
        try:
            super(SACTester, self).assert_test_info(*args, **kwds)
        finally:
            pass


if __name__ == "__main__":
    test_sac()
    # sact = SACTester("ShuffledAutocorrelation")
    # sact.audit = True
    # info = sact.run_test(plotflag=False)
    # print(info)
