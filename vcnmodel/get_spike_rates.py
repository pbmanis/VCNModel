"""
Compute mean input rates for each cell based on the revcorr data sets


This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2019 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
from pathlib import Path

import pickle

def get_rates():
    stim = 'Spont'
    fspont = f"GradeA_RCD_RCP_all_revcorrs_{stim:s}.pkl"
    if Path(fspont).is_file():
        print('is file')
    else:
        print('file not found')
    with open(fspont, 'rb') as fh:
        d = pickle.load(fh)
    print(d.keys())
    for c in d.keys():
        print(f"Cell: {c:6d}\n{'Source':<6s}  {'Interval (ms)':^14s} {'Rate (Hz)':^10s} {'N Spikes':8s}")
        print(f"{'Post':<6s} {d[c][0].mean_post_intervals:14.3f} {1000./d[c][0].mean_post_intervals:10.2f} {d[c][0].npost_spikes:8d}")
        for i in range(len(d[c][0].mean_pre_intervals)):
            print( f"Pre {i:2d} {d[c][0].mean_pre_intervals[i]:14.3f} {1000./d[c][0].mean_pre_intervals[i]:10.2f} {len(d[c][0].pre_st[i]):8d}")

if __name__ == '__main__':
    get_rates()