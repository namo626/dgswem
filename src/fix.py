#!/usr/bin/env python3

import re

with open('rhs_dg_hydro.f90') as f, open('rhs_dg_hydro.F90', 'w') as outfile:
    line = f.read()
    res = (re.findall('\n\\s+&', line))
    fixed = re.sub('\n\\s+&', ' &\n', line)

    fixed = re.sub('\n\\s+\$', ' &\n', fixed)
    outfile.write(fixed)
