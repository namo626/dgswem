#!/usr/bin/env python3
import pytest
from util import run_serial, run_parallel, TEST
import os

'''
Add tests by specifying as a 4-tuple in test_list:
1. serial/parallel run type
2. the name of the directory containing fort.14, fort.15, etc.
3. absolute tolerance (m)
4. relative tolerance
'''
test_list = [
    (run_serial, "quarter_annular", 0.05, 0.01),
    (run_parallel, "quarter_annular", 0.05, 0.01),
    (run_serial, "wetdry", 0.01, 1e-3),
    (run_parallel, "wetdry", 0.01, 1e-3),
    (run_serial, "mass_conservation", 1e-3, 1e-4),
]

@pytest.mark.parametrize("runner, name, atol, rtol", test_list)
def test_case(runner, name, atol, rtol):
    runner(os.path.join(TEST, name), atol=atol, rtol=rtol)


if __name__ == "__main__":
    pytest.main(["-vv"])
