import numpy as np
from typing import Union, List

def flatten_spike_array(spike_times:Union[np.ndarray, List],
                        time_window:tuple=(0., 1.0),
                        isi_flag=False):
    """Make a spike array "flat" (e.g, just a list of spike times across all repetitions)
    Optionally, do this to make an ISI distribution, which is accumulated across all spike
    train trials (and which only has ISIs that are within a trial)
    The returned spike times are clipped to the time_window, but retain their original
    times.

    Args:
        spike_times (Union[np.ndarray, List]): a 2D array of spike times, or a list of arrays
        time_window (float, optional): Time to use as the 0 marker. Defaults to 0.0.
        max_time (float, optional): Maximum time to include (before subtracting zero time). Defaults to 1.0.
        isi_flag (bool, optional): return ISIs instead of spike times. Defaults to False.

    Returns:
        _type_: _description_
    """
    spf = []
    for x in spike_times:
        x = np.array(x)
        x = x[(x >= time_window[0]) & (x < time_window[1])]
        if isi_flag:
            x = np.diff(x)
        spf.extend(x)
    spike_times_flat = np.array(spf, dtype=object).ravel()
    return spike_times_flat