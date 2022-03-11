from vcnmodel.util.user_tester import UserTester
import numpy as np
from elephant import spike_train_correlation as ESTC
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_generation import homogeneous_poisson_process
import quantities as pq 
# from vcnmodel.analyzers import spikestatistics as SPKS  # from Brian
from vcnmodel.analyzers import reverse_correlation as RC

"""
Test different cross-correlation methods (and compare)
"""
class regular_process(object):
    def __init__(self, start, interval, duration, offset=0.0):
        self.times = pq.s*(offset + np.arange(start, duration, interval))

        
def test_cross_correlation():
   CrossCorrelationTester(key="CrossCorrelation")


def compute_cc_data(method='coherent'):
    if method not in ['coherent', 'offset', 'uncorrelated']:
        raise ValueError(f"The methods for testing cc is not known: '{str(method):s}'")
    # st1 = homogeneous_poisson_process(rate=20.0 * pq.Hz, t_start=0.0 * pq.s, t_stop=50.0 * pq.s)
    st1 = regular_process(
            0.010, 1.0, 0.010
        )  # all identical
    if method == 'coherent':
        st2 = st1
    elif method == 'offset':
        st2 = regular_process(
            0.010, 1.0, 0.010, offset = 0.005,
        )  # al
    elif method == 'uncorrelated':
        # st2 = homogeneous_poisson_process(rate=20.0 * pq.Hz, t_start=0.0 * pq.s, t_stop=50.0 * pq.s)
        st2 = regular_process((0.1/np.pi), 1.0, (0.1/2.713))
    # st1 = np.linspace(0., 1, 100)*pq.s
    # st2 = st1 + 0.005*pq.s
    # standard analysis parameters:
    bw = 0.20*1e-3*pq.s
    width = 30.0*1e-3*pq.s
    # use simple spike correlator for analysis
    cc_simple = RC.reverse_correlation(st1.times, st2.times, binwidth=bw,
        corrwindow=[-width, width], )
    return cc_simple, st1, st2

class CrossCorrelationTester(UserTester):
    def init(self, key):
        UserTester.__init__(self, key)
    
    def assert_test_info(self, *args, **kwds):
        try:
            super(CrossCorrelationTester, self).assert_test_info(*args, **kwds)
        finally:
            pass
    
    def run_test(self):
        for m in ["coherent", "offset", "uncorrelated"]:
            self.ccs, self.st1, self.st2 = compute_cc_data(method=m)
            self.show_result()
        # if self.audit:
        #     self.show_result()

    def show_result(self):
        from matplotlib import pyplot as mpl
        width = 0.1
        f, ax = mpl.subplots(2, 1)

        ax[0].eventplot([self.st1.times, self.st2.times], linelengths=0.75)
        cc_len = len(self.ccs[0])
        t = np.linspace(-width, width, cc_len, endpoint=True)
        print("t: ", t.shape)
        print("ccs: ", cc_len)
        ax[1].plot(t, self.ccs[0], "r-")
        ax[1].plot(t[int(len(t)/2):], self.ccs[0][int(len(t)/2):], 'b-')
        # mpl.xlabel("time (" + str(cc_hist[0].times.units) + ")")
        mpl.ylabel("cross-correlation histogram")
        mpl.axis("tight")
        mpl.show()
