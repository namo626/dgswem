#!/usr/bin/env python3

import re

def format_F(fname):
    try:
        with open(fname) as f, open(fname[:-1] + 'F90','w') as outfile:
            line = f.read()
            fixed = re.sub('\n\\s+&', ' &\n', line)
            fixed = re.sub('\n\\s+\\$', ' &\n', fixed)
            fixed = re.sub('^c', '!', fixed, flags=re.IGNORECASE|re.MULTILINE)

            outfile.write(fixed)

    except FileNotFoundError:
        print("Error: input file %s not found." % (fname))


files = [
    'precipitation.F',
    'ocean_edge_hydro.F',
]

if __name__ == '__main__':
    for f in files:
        format_F(f)
