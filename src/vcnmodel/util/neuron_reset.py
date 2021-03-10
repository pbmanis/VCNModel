import numpy as np
import numpy.ma as ma # masked array
import re, sys, gc, collections

import neuron

"""
this is from cnmodel.util.pynrnutilities
Reset the neuron state. 
Only used in test modules (and then, only needed in those that run
some neuron code)

"""


_mechtype_cache = None
def all_mechanism_types():
    """Return a dictionary of all available mechanism types.
    
    Each dictionary key is the name of a mechanism and each value is
    another dictionary containing information about the mechanism::
    
        mechanism_types = {
            'mech_name1': {
                'point_process': bool,
                'artificial_cell': bool,
                'netcon_target': bool,
                'has_netevent': bool,
                'internal_type': int,
                'globals': {name:size, ...},
                'parameters': {name:size, ...},
                'assigned': {name:size, ...},
                'state': {name:size, ...},
            },
            'mech_name2': {...},
            'mech_name3': {...},
            ...
        }
    
    * point_process: False for distributed mechanisms, True for point 
        processes and artificial cells.
    * artificial_cell: True for artificial cells, False otherwise
    * netcon_target: True if the mechanism can receive NetCon events
    * has_netevent: True if the mechanism can emit NetCon events
    * internal_type: Integer specifying the NEURON internal type index of 
        the mechanism
    * globals: dict of the name and vector size of the mechanism's global
        variables
    * parameters: dict of the name and vector size of the mechanism's 
        parameter variables
    * assigned: dict of the name and vector size of the mechanism's 
        assigned variables
    * state: dict of the name and vector size of the mechanism's state
        variables


    Note: The returned data structure is cached; do not modify it.
        
    For more information on global, parameter, assigned, and state 
    variables see:
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/nmodl/nmodl.html
    """
    global _mechtype_cache
    if _mechtype_cache is None:
        _mechtype_cache = collections.OrderedDict()
        mname = neuron.h.ref('')
        # Iterate over two mechanism types (distributed, point/artificial)
        for i in [0, 1]:
            mt = neuron.h.MechanismType(i)
            nmech = int(mt.count())
            # Iterate over all mechanisms of this type
            for j in range(nmech):
                mt.select(j)
                mt.selected(mname)
                
                # General mechanism properties
                name = mname[0]  # convert hoc string ptr to python str
                
                desc = {
                    'point_process': bool(i),
                    'netcon_target': bool(mt.is_netcon_target(j)),
                    'has_netevent': bool(mt.has_net_event(j)),
                    'artificial_cell': bool(mt.is_artificial(j)),
                    'internal_type': int(mt.internal_type()),
                }
                
                # Collect information about 4 different types of variables
                for k,ptype in [(-1, 'globals'), (1, 'parameters'), 
                                (2, 'assigned'), (3, 'state')]:
                    desc[ptype] = {} # collections.OrderedDict()
                    ms = neuron.h.MechanismStandard(name, k)
                    for l in range(int(ms.count())):
                        psize = ms.name(mname, l)
                        pname = mname[0]  # parameter name
                        desc[ptype][pname] = int(psize)
                
                # Assemble everything in one place
                _mechtype_cache[name] = desc
            
    return _mechtype_cache

def reset(raiseError=True):
    """Introspect the NEURON kernel to verify that no objects are left over
    from previous simulation runs.
    """
    # Release objects held by an internal buffer
    # See https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=3221
    neuron.h.Vector().size()    
    
    # Make sure nothing is hanging around in an old exception or because of
    # reference cycles 

    # sys.exc_clear()
    gc.collect(2)
    neuron.h.Vector().size()
    numsec = 0

    remaining = []
    n = len(list(neuron.h.allsec()))

    if n > 0:
        remaining.append((n, 'Section'))
        
    n = len(neuron.h.List('NetCon'))
    if n > 0:
        remaining.append((n, 'NetCon'))
    
    # No point processes or artificial cells left
    for name, typ in all_mechanism_types().items():
        if typ['artificial_cell'] or typ['point_process']:
            n = len(neuron.h.List(name))
            if n > 0:
                remaining.append((n, name))
    
    if len(remaining) > 0 and raiseError:  # note that not raising the error leads to memory leak
        msg = ("Cannot reset--old objects have not been cleared: %s" %
               ', '.join(['%d %s' % rem for rem in remaining]))
        raise RuntimeError(msg)
