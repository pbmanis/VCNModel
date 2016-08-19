__author__ = 'pbmanis'


import multiprocessing
import time
import model_run
import itertools
import pickle
import numpy as np

celltype = 'Stellate'
model_type = 'XM13'

start_time = time.time()

def run_model_Star(*args):
    """
    This parses out the parameter exploration, for parallel processing
    """
    thisModel = model_run.ModelRun()
    thisModel.set_celltype(celltype)
    thisModel.set_model_type(model_type)
    thisModel.set_starttime(start_time)
    return thisModel.run_model(*args)

# define parameter space to explore
#
#
if celltype == 'Stellate':
    space = {'ihvcn': np.linspace(0.2, 2.0, 19), 'kht': np.linspace(0.5, 2.0, 19)}
elif celltype == 'Bushy':
    space = {'ihvcn': np.linspace(0.2, 2, 5), 'klt': np.linspace(0.2, 2, 5)}
else:
    print 'What cell type? %s' % (celltype)
    exit()

spvalues = [x for x in apply(itertools.product, space.values())]
spaceMap = [dict(zip(space.keys(), p)) for p in spvalues]
for idnum, x in enumerate(spaceMap): # add an ID as a key
    x['id'] = idnum

test = False  # false for parallel runs, true for serial testing
######
#  non parallel run with no save: testing for errors
if test is True:
    for s in spaceMap:
        run_model_Star(s)

    runtime = time.time() - start_time
    print 'run time: ', runtime, "seconds"
    exit()
######

PROCESSES = len(spaceMap)
TASKS = [s for s in spaceMap]
print "Number of tasks: ", len(TASKS)

ncpu = multiprocessing.cpu_count()
print 'ncpu: ', ncpu
pool = multiprocessing.Pool(ncpu)

if len(TASKS) < ncpu:
    chunksize = len(TASKS)
else:
    chunksize = len(TASKS) // ncpu
print 'chunksize: ', chunksize

#print 'tasks: ', TASKS


imap_it = pool.imap(run_model_Star, TASKS, chunksize)

n = 0
result = {}
#try:
for x in imap_it:
    if x is not None:
        result[n] = x
        #print 'Result %d:\n' % (n), result[n]
        n = n + 1
    else:
        print '>'*80
        print '   parallel_model_run: multiprocessing result error: X is NONE?'
#except:
#    print '$'*80
#    print 'try FAILED at n = %d' % n

pool.close()
pool.join()

# print '*'*80
# print 'result: ', result
# print '*'*80

big_result = {'space': space, 'spacemap': spaceMap, 'result': result}
fn = celltype + '_' + model_type + '_'
for k in space.keys():
    if k == 'id':
        continue
    fn = fn + k + '_' # add variables for run
#dtime = time.strftime("%y.%m.%d-%H.%M.%S")
dtime = time.strftime("%y.%m.%d-%H.%M.%S", time.localtime(start_time))  # use run starttime
fn += dtime + '.p'
f = open(fn, 'wb')
pickle.dump(big_result, f)
f.close()

runtime = time.time() - start_time
print 'run time: ', runtime, "seconds"
print 'file: ', fn
