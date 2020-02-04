"""
Parse an swc
"""

import re
import argparse
from pathlib import Path


re_section = re.compile('\s*(sections\[)([0-9]*)\]\s*{')
re.compile('\s*(pt3dadd\()([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)\,\s([-+]?[0-9]*\.?[0-9]+)')
re_section = re.compile('\s*(sections\[)([0-9]*)\]\s*{')
re_access = re.compile('\s*(access)\s*(sections\[)([0-9]*)\]\s*')
re_append = re.compile('\s*([a-z]*)(.append\(\))')
re_connect = re.compile('\s*(connect)\s*(sections\[)([0-9]*)\](\([0-9]*\)),\s*(sections\[)([0-9]*)\](\([0-9]*\))')

re_seg = re.compile('(seg\=)([\d]*)') # '([d+])$')
re_endsec = re.compile('^}')


def make_map(fn):
    dout = ''
    in_section = False
    secstr = ''
    with open(fn, 'r') as fh: 
        for cnt, line in enumerate(fh):  # read the input file line by line
            line = line.rstrip().lstrip()
            s = re_section.match(line)
            if s is not None:
                secno = s.groups()[1]
                secstr = (f'section[{secno:s}]: ')
                in_section = True
                swcs = []
                continue
            if in_section:
                if re_endsec.match(line):
                    in_section = False
                    for swi in swcs:
                        secstr +=(f"{swi:s}, ")
                    print(secstr)
                    dout += secstr + '\n'
                    secstr = [] # reset
                    continue
                swcindex = re_seg.search(line)
                if swcindex is not None:
                    swci = swcindex.groups()[1]
                    swcs.append(swci)
    fout = Path(Path(fn).name).with_suffix('.segmap')
    fout.write_text(dout)


def main():
    # fn = 'VCN_Cells/VCN_c09/Morphology/VCN_c09_FullCell_edited.hoc'
    parser = argparse.ArgumentParser(description='hoc_swc_sectionmapping',
                    argument_default=argparse.SUPPRESS,
                    fromfile_prefix_chars='@')
    parser.add_argument(dest='input_file', action='store',
                   default=None,
                   help='Select the hoc file to map (no default)')
    args = parser.parse_args()
    make_map(args.input_file)

if __name__ == '__main__':
    main()
