__author__ = 'pbmanis'

import pyqtgraph as pg

class Render():
    def __init__(self, hf):
        self.hf = hf

    def get_hoc_file(self, infile):
        if self.hf.file_loaded is False:
            exit()
        if not HAVE_PG:
            return()
        self.section_list = self.hf.get_section_prefixes()
        self.hf.sec_groups.keys()
        if len(self.hf.sec_groups) > 1: # multiple names, so assign colors to structure type
            self.section_colors = {}
            for i, s in enumerate(self.hf.sec_groups.keys()):
                self.section_colors[s] = self.hg.get_color_map(i)
        else: # single section name, assign colors to SectionList types:
            self.section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
                'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}

        (v, e) = self.hf.get_geometry()
        self.clist = []

        for si in self.hf.sections: # self.section_list[s]:
            self.hf.h('access %s' % si)
            sr = self.hf.h.SectionRef()
            n1 = self.hf.h.cas().name()
            if sr.has_parent() == 1:
                x=sr.parent
                n2 = x.name()
                self.clist.append([n1, n2])
            else:
                self.clist.append([n1, None])

    def render(self, mech):
        pg.mkQApp()
        render = HocViewer(self.hf)

        type = 'line'
        if type == 'line':
            line = render.draw_graph()
            line.set_group_colors(self.section_colors, alpha=0.35)
        if type == 'surface':
            surface = render.draw_surface(resolution = 1.0)
            surface.set_group_colors(self.section_colors, alpha=0.35)
        elif type == 'cylinder':
            cylinder=render.draw_cylinders()
            cylinder.set_group_colors(self.section_colors, alpha=0.35)
        elif type == 'volume':
            volume = render.draw_volume(resolution = 1.0, max_size=1e9)
            #volume.set_group_colors(section_colors, alpha=0.35)

    #        render.hr.read_hoc_section_lists(self.section_colors.keys())
        #surface = render.draw_surface()
        # cylinder = render.draw_cylinders()
        # cylinder.set_group_colors(self.section_colors, alpha=0.35, mechanism=mech)
        # line = render.draw_graph()
       # surface.set_group_colors(self.section_colors, alpha=0.35)
       # surface.set_group_colors(self.section_colors, alpha=0.35, mechanism=mech)
