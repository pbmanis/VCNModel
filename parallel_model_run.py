__author__ = 'pbmanis'


import multiprocessing
import time
import model_run
import itertools
import pickle
import numpy as np

celltype = 'Stellate'
modeltype = 'XM13'

def run_model_Star(*args):
    """
    This parses out the parameter exploration, for parallel processing
    """
    thisModel = model_run.ModelRun()
    thisModel.set_celltype(celltype)
    thisModel.set_modeltype(modeltype)
    return thisModel.runModel(*args)

start_time = time.time()
# define parameter space to explore
# use logspace with odd number of values so that middle is always 1x
#
space = {'ihvcn': np.logspace(0.1, 10, 5), 'kht': np.logspace(0.1, 10, 5)}
spvalues = [x for x in apply(itertools.product, space.values())]
spaceMap = [dict(zip(space.keys(), p)) for p in spvalues]


# non parallel run with no save: testing for errors
for s in spaceMap:
    run_model_Star(s)

exit()
#

PROCESSES = len(spaceMap)
pool = multiprocessing.Pool(PROCESSES)

TASKS = [vals for vals in spaceMap]
print "Number of tasks: ", len(TASKS)

imap_it = pool.imap(run_model_Star, TASKS)

n = 0
result = {}
#try:
for x in imap_it:
    if x is not None:
        result[n] = x
        print 'Result %d:\n' %(n) , result[n]
        n = n + 1
    else:
        print '>'*80
        print '   parallel_model_run: multiprocessing result error: X is NONE?'
#except:
#    print '$'*80
#    print 'try FAILED at n = %d' % n

pool.close()
pool.join()

print '*'*80
print 'result: ', result
print '*'*80

big_result = {'space': space, 'spacemap': spaceMap, 'result': result}
fn = celltype + '_' + modeltype + '_'
fn = fn + ''.join([k+'_' for k in space.keys()])  # add variables for run
dtime = time.strftime("%y.%m.%d-%H.%M.%S")
fn += dtime + '.p'
f = open(fn, 'wb')
pickle.dump(big_result, f)
f.close()

runtime = time.time() - start_time
print 'run time: ', runtime, "seconds"

