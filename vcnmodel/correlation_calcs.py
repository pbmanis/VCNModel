import numpy as np
import scipy
from typing import Dict, List, Union, Tuple
from pylibrary.tools.cprint import cprint as CP

"""
This code was taken from McTavish's code.

Main change is that we specify the CCF window as tstart--tend. rather than a window around
0.

"""
"""
Spike train utility methods.

Utility methods for reading spike train files as well as to analyze spike train 
matrices and vectors.

.. seealso::
    :cmv colass:`~neuronpy.graphics.spikeplot.SpikePlot` for visualization of spike 
    trains.

AUTHORS:

- THOMAS MCTAVISH (2010-03-01): initial version
- THOMAS MCTAVISH (2011-05-03): Additions mostly for synchrony analysis and 
    filtering. Refactored as ``spiketrain`` instead of ``spiketrainutil``.
"""
# While this software is under the permissive MIT License,
# (http://www.opensource.org/licenses/mit-license.php)
# We ask that you cite the neuronpy package (or tools used in this package)
# in any publications and contact the author with your referenced publication.
#
# Format:
# McTavish, T.S. NeuronPy library, version 0.1,
# http://bitbucket.org/tommctavish/neuronpy
#
# Copyright (c) 2010 Thomas S. McTavish
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

def get_sync_masks(train_a:Union[np.ndarray, List],
                   train_b:Union[np.ndarray, List], 
                   window:Union[List, np.ndarray, Tuple]=[-50., 50.], ) -> Tuple:
    """
    For two spike trains, return the mask of those spikes that
    are within some time window of co-occurrence in the other train.
    
    :param train_a: A list of spike times.
    
    :param train_b: Another list of spike times.
    
    :param window: Time window +/- about a given spike in one train to look
        for a co-occuring spike in the other train.
        
    :return: Two vectors of ``len(train_a)`` and ``len(train_b)`` where a
        zero indicates that a spike does not co-occur in the other train, and
        1 indicates that a spike co-occurs in the other train within 
        ``window`` time.
    """
    idx_a = 0
    idx_b = 0
    
    mask_a = np.zeros_like(train_a)
    mask_b = np.zeros_like(train_b)
    
    len_a = len(train_a)
    len_b = len(train_b)
    
    while idx_a < len_a and idx_b < len_b:
        val_a = train_a[idx_a]
        val_b = train_b[idx_b]

        diff = val_a - val_b
        if window[0] <= diff <= window[1]:
            mask_a[idx_a] = 1
            mask_b[idx_b] = 1
        
        if val_a == val_b:
            idx_a += 1
            idx_b += 1
        else:
            if val_a < val_b:
                idx_a += 1
            else:
                idx_b += 1
    
    return mask_a, mask_b

def get_sync_traits(train_a:Union[np.ndarray, List],
                   train_b:Union[np.ndarray, List], 
                   window:Union[List, np.ndarray, Tuple]=[-50., 50.],):
    """
    For two spike trains, get their masks where they have spikes that occur
    within ``window`` time of each other and the ratio of correlated vs. 
    total spikes in both trains.
    
    :param train_a: List of spike times in one train.
    
    :param train_b: List of spike times in another train.
    
    :param window: Time window to search for correlated inputs.
    
    :return: mask_a, mask_b, ratio -- the correlated masks and the ratio of
        correlated vs. total spikes in both trains.
    """
    mask_a, mask_b = get_sync_masks(train_a, train_b, window)
    len_a = len(train_a)
    len_b = len(train_b)
    num_coincident = np.sum(mask_a)
    ratio = 2. * float(num_coincident) / float(len_a + len_b)
    
    return num_coincident, mask_a, mask_b, ratio
        
def closest_timing(reference_train:Union[np.ndarray, List],
                   comparing_train:Union[np.ndarray, List], 
                   window:Union[List, np.ndarray, Tuple]=[-50., 50.],):
    """
    For each spike in the reference train, determine the closest spike time
    in the other train at least within some window.
    
    :param reference_train: List of spike times of the reference train.
    :param comparing_train: List of spike times in the comparing train.
    :param window: Time window in ms to search.
    
    :return: A dict where the keys are indices of the reference train and the
            time difference of the closest spike in the comparing train is the
            value.
    """
    time_dict = {}
    try:
        ita = reference_train.__iter__()
        itb = comparing_train.__iter__()

        val_a = next(ita)
        val_b = next(itb)
        idx_a = 0
        idx_b = 0
        while(True):
            diff = val_b - val_a
            print('diff: ', diff, window)
            if window[0] > diff  or window[1] < diff:
                print("  not  IN WINDOW")
                if val_a > val_b:
                    val_b = next(itb)
                    idx_b += 1
                else:
                    val_a = next(ita)
                    idx_a += 1
            else:
            # ***************************
                print('IN WINDOW')
                if idx_a in time_dict:
                    if diff < np.abs(time_dict[idx_a]):
                        time_dict[idx_a] = val_a - val_b
                else:
                    time_dict[idx_a] = val_a - val_b
                    
                if val_a == val_b:
                    val_a = next(ita)
                    idx_a += 1
                    val_b = next(itb)
                    idx_b += 1
                else:
                    if val_a > val_b:
                        val_a = next(ita)
                        idx_a += 1
                    else:
                        val_b = next(itb)
                        idx_b += 1
    except StopIteration:
        return time_dict
    
    return time_dict
    

def coincidence_factor(ref:Union[List, np.ndarray], 
            comp:Union[List, np.ndarray], 
            window:Union[List, np.ndarray, Tuple]=[-5., 1.], 
            isi:Union[float, None]=None) -> float:
    r"""
    The coincidence factor :math:`\Gamma` between two spike trains is defined as

    .. math::
        
       \Gamma = \frac{N_\mathrm{coinc}- E \left( N_\mathrm{coinc} \right)}
       {\frac{1}{2}\left(N_\mathrm{ref}+N_\mathrm{comp}\right) - 
       E \left( N_\mathrm{coinc} \right)}

    where :math:`N_{\mathrm{ref}}` are the number of spikes in the reference train,
    :math:`N_{\mathrm{comp}}` is the number of spikes in the comparing train, 
    :math:`N_{\mathrm{coinc}}` is the number of coincident spikes within a time window 
    :math:`\Delta`, :math:`E \left( N_\mathrm{coinc} \right) = 2 v \Delta N_{\mathrm{ref}}` 
    is the expected number of coincident spikes that would be given by chance 
    if the spikes in the comparing train were generated by a homogeneous 
    Poisson process with its rate :math:`v`. This correlation measure has the range 
    [-1, 1] where 1 is perfectly correlated, 0 is not correlated, and -1 is 
    perfectly anti-correlated.
    
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are synchronized.
        This has a default value of 5.
        
    :param isi: If supplied, this is the isi of the comparing train. Otherwise,
        the rate of the train is computed by taking the last spike minus the
        first spike and dividing by the number of spikes in between.
    
    :return: Coincidence factor
    """
    num_coincident, mask_a, mask_b, ratio = get_sync_traits(ref, comp, window)
    len_ref = len(ref)
    len_comp = len(comp)
    total_spikes = len_ref + len_comp
    if isi is None:
        v = (len_comp - 1) / (comp[-1] - comp[0])
    else:
        v = 1.0 / isi
    wdur = window[1]-window[0]
    expected_coincidences = 2 * v * wdur * len_ref
    return (
        (num_coincident - expected_coincidences)
        * 2
        / (total_spikes - (2 * expected_coincidences))
    )


def coincidence_factor_phase(ref:Union[List, np.ndarray], 
            comp:Union[List, np.ndarray], 
            window:Union[List, np.ndarray, Tuple]=[-5., 1.], 
            num_intervals:int=13, isi:Union[float, None]=None) -> Union[List, np.ndarray]:
    """
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the comparing train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are synchronized.
        This has a default value of 5.
        
    :param num_intervals: Number of iterations to perform, sliding
        ``comp`` from [-1/2, 1/2] the median of ``ref``'s
        interspike interval. This should be an odd number to ensure a 
        precise sample about 0 delay.
    
    :param isi: If supplied, this is the isi of the comparing train. Otherwise,
        the rate of the train is computed by taking the last spike minus the
        first spike and dividing by the number of spikes in between.
        
    :return: A vector of length ``num_intervals`` that corresponds to 
        coincidence factor values from a shift of -isi/2 to isi/2.
    """
    phi_vec = np.zeros(num_intervals)
    idx = 0
    if isi is None or not isinstance(isi, float):
        isi = get_mean_isi(comp)
    shift = isi / num_intervals / 2.0
    for i in np.linspace(-(isi / 2.0) + shift, (isi / 2.0) - shift, num_intervals):
        vec = np.add(comp, i)
        phi_vec[idx] = coincidence_factor(ref, vec, window, isi)
        idx += 1

    return phi_vec


def coincident_spikes(ref:Union[List, np.ndarray], 
            comp:Union[List, np.ndarray], 
            window:Union[List, np.ndarray, Tuple]=[-5., 1.],
            normalize:bool=False, 
            compfreq:Union[None, float]=None) -> Union[List, np.ndarray]:
    """
    Get the fraction of coincident spikes between two trains. Coincidence is
    defined as a spike co-occurring within some time window.
    
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the comparing train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are coincident.
        This has a default value of 5 ms.
    
    :param normalize: If ``True``, values are normalized by the rate expected by
        chance, which is defined as ``2 * frequency(comp) * window * len(ref)``.
        Also, a tuple is returned as 
        (normalized_coincidences, coincidences, expected_coincidences, total_spikes)
        where total_spikes is the length of both trains.

    :param compfreq: Frequency in Hz of comparing train. If None, the mean frequency
        of the comparing train is calcluated and used. This is only used when 
        *normalize* is ``True``.
        
    :return: A vector of length ``1 + 2*timewindow/dt``. If normalize is ``True``,
        return a tuple as 
        (normalized_coincidences, coincidences, expected_coincidences, total_spikes).
    """
    coincidences, mask_a, mask_b, ratio = get_sync_traits(ref, comp, window)
    if normalize == False:
        return coincidences

    len_ref = len(ref)
    len_comp = len(comp)
    total_spikes = len_ref + len_comp
    if compfreq is None:
        compfreq = (len_comp - 1) / (comp[-1] - comp[0])
    wdur = window[1]-window[0]
    expected_coincidences = 2 * compfreq * wdur * len_ref
    return (
        float(coincidences) / float(expected_coincidences),
        coincidences,
        expected_coincidences,
        total_spikes,
    )


def coincident_spikes_correlogram(
                ref:Union[List, np.ndarray]=None, 
                comp:Union[List, np.ndarray]=None, 
                window:Union[List, np.ndarray, Tuple]=[-5., 1.],
                binwidth:float=0.1,
                timewindow:Union[List, Tuple]=[-100., 100.], 
                normalize:bool=False
                ) -> Union[np.ndarray, List]:
    """
    Get the correlogram of coincident spikes between two trains. Coincidence is
    defined as a spike co-occurring within some time window. The number of 
    coincidences is returned as a vector of len(1 + (2*timewindow)). This means that
    the middle index is the value at lag zero ms.
    
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the comparing train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are synchronized.
        This has a default value of 5.
        
    :param dt: The binwidth.
    
    :param timewindow: Correlogram range between [-timewindow, timewindow].
        
    :param normalize: If True, values are normalized by the rate expected by
        chance at each lag, where chance is defined as 
        ``2 * frequency(comp) * window * len(ref)``.
        
    :return: A vector of length ``1 + 2*timewindow/dt``.
    """
    if ref == None or comp == None:
        raise ValueError('coincident_spikes_correlation: reference and comparator must be defined')
    phi_vec = np.zeros(int(1 + ((timewindow[1]-timewindow[0]) / float(binwidth))))
    for idx in range(len(phi_vec)):
        i = -(timewindow[1]-timewindow[0]) + (idx * binwidth)
        vec = np.add(comp, i)
        result = coincident_spikes(ref, vec, window, normalize=normalize)
        if normalize:
            phi_vec[idx] = result[0]
        else:
            phi_vec[idx] = result

    return phi_vec


def revcorr_spikes(ref:Union[List, np.ndarray]=None, 
                comp:Union[List, np.ndarray]=None, 
                window:Union[List, np.ndarray, Tuple]=[-5., 1.],
                binwidth:float=0.1,
                timewindow:Union[List, Tuple]=[-100., 100.], 
                ):
    if ref is None or comp is None:
        raise ValueError("coincident_spikes_correlation: reference and comparator must be defined")
    refa = [a for a in ref if (timewindow[0] <= a <= timewindow[1])]
    refb = [b for b in comp if (timewindow[0] <= b <= timewindow[1])]
    td = closest_timing(refa, refb, window)
    # print('refa: ', refa)
    # print('refb: ', refb)
    # print('td: ', td)
    return td
    

def test_revcorr():
    tb = np.arange(0, 100., 0.1)
    st1 = []
    st2 = []
    print(len(tb))
    for i in range(0, len(tb), 10):
        st1.append(tb[i])
        j = i
        if j < len(tb):
            st2.append(tb[j])
    print(st1)
    print(st2)
    rcd = coincident_spikes_correlogram(st1, st2, window=[-10., 10.],
        timewindow = [0., np.max(tb)])
    print(rcd)
    rc = rcd # list(rcd.values())
    import matplotlib.pyplot as mpl

    rci = [int(r/0.1) for r in rc]
    rcdata = np.zeros_like(tb)
    # rcdata[rci] = rc[rci]
    mpl.plot(tb, rcd)
    mpl.show()



if __name__ == '__main__':
    test_revcorr()
    
    