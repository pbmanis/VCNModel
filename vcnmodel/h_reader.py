from __future__ import print_function
__author__ = 'pbmanis'


try:
    from neuron import h
    import neuron
    HAVE_NEURON=True
except:
    HAVE_NEURON=False

import collections
import numpy as np

import os
import re

import os.path
class HocReader(object):
    """
    Provides useful methods for reading hoc structures.

    Input:
        hoc: a hoc object or a "xxx.hoc" file name.
    """
    def __init__(self, hoc):
        self.file_loaded = False
        if isinstance(hoc, basestring):
            fullfile = os.path.join(os.getcwd(), hoc)
            if not os.path.isfile(fullfile):
                raise Exception("File not found: %s" % (fullfile))
            success = neuron.h.load_file(1, fullfile)
            if success == 0: # indicates failure to read the file
                raise NameError("Found file, but NEURON load failed: %s" % (fullfile))
            self.file_loaded = True
            self.h = h # save a copy of the hoc object itself.
        else:
            self.h = hoc # just use the passed argument
            self.file_loaded = True

        # geometry containers
        self.edges = None
        self.vertexes = None

        # all sections in the hoc  {sec_name: hoc Section}
        self.sections = collections.OrderedDict()
        # {sec_name: index} indicates index into self.sections.values() where
        # Section can be found.
        self.sec_index = {}
        # {sec_name: [mechanism, ...]}
        self.mechanisms = collections.OrderedDict()
        # {sec_group_name: set(sec_name, ...)}
        self.sec_groups = {}

        # populate self.sections and self.mechanisms
        self._read_section_info()

        # auto-generate section groups based on either hoc section lists, or
        # on section name prefixes.
        sec_lists = self.get_section_lists()
        sec_prefixes = self.get_section_prefixes()


        # Add groupings by section list if possible:
        if len(sec_lists) > 1:
            self.add_groups_by_section_list(sec_lists)

        # Otherwise, try section prefixes
        elif len(sec_prefixes) > 1:
            for group, sections in sec_prefixes.items():
                self.add_section_group(group, sections)


    def update(self):
        """
        Update information on sections after external changes
        """
        self._read_section_info()


    def get_section(self, sec_name):
        """
        Return the hoc Section object with the given name.
        """
        try:
            return self.sections[sec_name]
        except KeyError:
            raise KeyError("No section named '%s'" % sec_name)


    def get_section_prefixes(self):
        """
        Go through all the sections and generate a dictionary mapping their
        name prefixes to the list of sections with that prefix.

        For example, with sections names axon[0], axon[1], ais[0], and soma[0],
        we would generate the following structure:

            {'axon': ['axon[0]', 'axon[1]'],
             'ais':  ['ais[0]'],
             'soma': ['soma[0]']}
        """
        prefixes = {}
        regex = re.compile('(?P<prefix>\w+)\[(\d*)\]')
        for sec_name in self.sections:
            g = regex.match(sec_name)
            if g is None:
                continue
            prefix = g.group('prefix')
            prefixes.setdefault(prefix, []).append(sec_name)
        return prefixes


    def get_mechanisms(self, section):
        """
        Get a set of all of the mechanisms inserted into a given section
        Input: section (hoc object)
        Returns: set of mechanism names
        Side-effects: None
        """
        return self.mechanisms[section]


    def get_density(self, section, mechanism):
        """
        Get density mechanism that may be found the section.
        mechanism is a list ['name', 'gbarname']. This is needed because
        some mechanisms do not adhere to any convention and may have a different
        kind of 'gbarname' than 'gbar<name>_mechname'
        returns the average of the conductance density, as that may range across different
        values in a section (e.g., can vary by segments)
        Input: section (hoc object)
                mechanism mechanism is a list ['name', 'gbarname'].
        Output:
            mean conductance inserted into the section across segments
        Side-effects:
            None
        """

        gmech = []
        for seg in section:
            try:
                x =  getattr(seg,  mechanism[0])
                mecbar = '%s_%s' % (mechanism[1], mechanism[0])
                if mecbar in dir(x):
                    gmech.append(getattr(x, mechanism[1]))
                else:
                    print('hoc_reader:get_density did not find the mechanism in dir x')
            except NameError:
                return(0.)
            except:
                print('hoc_reader:get_density failed to evaluate the mechanisms... ')
                raise


        print(gmech)
        if len(gmech) == 0:
            gmech = 0.
        return np.mean(gmech)


    def get_sec_info(self, section):
        """
        Get the info of the given section
        modified from: neuronvisio
        """
        info = "<b>Section Name:</b> %s<br/>" %section.name()
        info += "<b>Length [um]:</b> %f<br/>" % section.L
        info += "<b>Diameter [um]:</b> %f<br/>" % section.diam
        info += "<b>Membrane Capacitance:</b> %f<br/>" % section.cm
        info += "<b>Axial Resistance :</b> %f<br/>" % section.Ra
        info += "<b>Number of Segments:</b> %f<br/>" % section.nseg
        mechs = []
        for seg in section:
            for mech in seg:
                mechs.append(mech.name())
        mechs = set(mechs) # Excluding the repeating ones

        mech_info = "<b>Mechanisms in the section</b><ul>"
        for mech_name in mechs:
            s = "<li> %s </li>" % mech_name
            mech_info += s
        mech_info += "</ul>"
        info += mech_info
        return info


    def _read_section_info(self):
        # Collect list of all sections and their mechanism names.
        self.sec_index=collections.OrderedDict()
        self.sections = collections.OrderedDict()
        self.mechanisms = collections.OrderedDict()
        for i, sec in enumerate(self.h.allsec()):
            self.sections[sec.name()] = sec
            self.sec_index[sec.name()] = i
            mechs = set()
            for seg in sec:
                for mech in seg:
                    mechs.add(mech.name())
            self.mechanisms[sec.name()] = mechs


    def hoc_namespace(self):
        """
        Return a dict of the HOC namespace {'variable_name': hoc_object}.
        NOTE: this method requires NEURON >= 7.3
        """
        names = {}
        for hvar in dir(self.h): # look through the whole list, no other way
            try:
                # some variables can't be pointed to...
                if hvar in ['nseg', 'diam_changed', 'nrn_shape_changed_',
                            'secondorder', 'stoprun']:
                    continue
                u = getattr(self.h, hvar)
                names[hvar] = u
            except:
                continue
        return names


    def find_hoc_hname(self, regex):
        """
        Return a list of the names of HOC objects whose *hname* matches regex.
        """
        objs = []
        ns = self.hoc_namespace()
        for n, v in ns.items():
            try:
                hname = v.hname()
                if re.match(regex, hname):
                    objs.append(n)
            except:
                continue
        return objs


            # m = sid.match(hname)
            # sections=[]
            # if m is not None:
            #     for v in getattr(self.h, hvar):
            #         sections.append(v)
            #     self.add_section_group(hvar, sections)


    def add_section_group(self, name, sections, overwrite=False):
        """
        Declare a grouping of sections (or section names). Sections may be
        grouped by any arbitrary criteria (cell, anatomical type, etc).

        Input:
            name: string name of the section group
            sections: list of section names or hoc Section objects.

        """
        if name in self.sec_groups and not overwrite:
            raise Exception("Group name %s is already used (use overwrite=True)." % name)

        group = set()
        for sec in sections:
            if not isinstance(sec, basestring):
                sec = sec.name()
            group.add(sec)
        self.sec_groups[name] = group


    def get_section_group(self, name):
        """
        Return the set of section names in the group *name*.
        """
        return self.sec_groups[name]


    def get_section_lists(self):
        """
        Search through all of the hoc variables to find those that are "SectionLists"
        """
        return self.find_hoc_hname(regex=r'SectionList\[')
        #ns = self.hoc_namespace()
        #return [name for name in ns if ns[name].hname().startswith('SectionList[')]


    def add_groups_by_section_list(self, names):
        """
        Add a new section groups from the hoc variables indicated in *names*.

        Input:
            names: list of variable names. Each name must refer to a list of
                   Sections in hoc. If a dict is supplied instead, then it
                   maps {hoc_list_name: section_group_name}.
        Side effects (modifies):
           calls add_section_group
        returns: Nothing.
        """
        # if a list is supplied, then the names of groups to create are
        # exactly the same as the names of hoc lists.
        if not isinstance(names, dict):
            names = {name:name for name in names}
        for hoc_name, group_name in names.items():
            var = getattr(self.h, hoc_name)
            self.add_section_group(group_name, list(var))


    def get_geometry(self):
        """
        modified from:neuronvisio
        Generate structures that describe the geometry of the sections and their segments (all segments are returned)
        Inputs: None
        Outputs: vertexes: record array containing {pos: (x,y,z), dia, sec_id}
                           for each segment.
                 edges:    array of pairs indicating the indexes of connected
                           vertexes.
        Side effects:
            modifies vertexes and edges.
        """

        # return cached geometry if this method has already run.
        if self.vertexes is not None:
            return self.vertexes, self.edges

        self.h.define_shape()

        # map segments (lines) to the section that contains them
        self.segment_to_section = {}

        vertexes = []
        connections = []

        for secid, sec in enumerate(self.sections.values()):
            x_sec, y_sec, z_sec, d_sec = self.retrieve_coordinate(sec)

            for i,xi in enumerate(x_sec):
                vertexes.append(((x_sec[i], y_sec[i], z_sec[i]), d_sec[i], secid))
                indx_geom_seg = len(vertexes) - 1
                if len(vertexes) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])

        self.edges = np.array(connections)
        self.vertexes = np.empty(len(vertexes), dtype=[
            ('pos', float, 3),
            ('dia', float),
            ('sec_index', int)])
        self.vertexes[:] = vertexes
        return self.vertexes, self.edges


    def retrieve_coordinate(self, sec):
        """Retrieve the coordinates of the section avoiding duplicates"""

        sec.push()
        x, y, z, d = [],[],[],[]

        tot_points = 0
        connect_next = False
        for i in range(int(self.h.n3d())):
            present = False
            x_i = self.h.x3d(i)
            y_i = self.h.y3d(i)
            z_i = self.h.z3d(i)
            d_i = self.h.diam3d(i)
            # Avoiding duplicates in the sec
            if x_i in x:
                ind = len(x) - 1 - x[::-1].index(x_i) # Getting the index of last value
                if y_i == y[ind]:
                    if z_i == z[ind]:
                        present = True

            if not present:
                k =(x_i, y_i, z_i)
                x.append(x_i)
                y.append(y_i)
                z.append(z_i)
                d.append(d_i)
        self.h.pop_section()
        return (np.array(x),np.array(y),np.array(z),np.array(d))

