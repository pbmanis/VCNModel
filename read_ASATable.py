"""
Read the ASA table (by size), as exported into a csv file
"""

import pandas as pd
import numpy as np

fn = 'VCN_ASA.csv'

df = pd.read_csv(fn)
print df['7'][0]

for k in df.keys():
    n = 0 # count how many
    while not np.isnan(df[k][n]):
        if n == 0:
            print 'VCN_c%02d = [ ' % int(k)
        print "    [(%.2f), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ]," % df[k][n]
        n += 1
    
    if n > 0:
        print '    ]\n'

