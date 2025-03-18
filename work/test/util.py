#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import subprocess
import pytest
import os
import git


repo = git.Repo('.', search_parent_directories=True)
ROOT = repo.working_tree_dir
WORK = os.path.join(ROOT, "work")
TEST = os.path.join(WORK, "test")

def mean_rel_err(xs, ys):
    ab = np.abs(xs - ys)
    mask = ab < 90000
    e = ab[mask]
    e = e / np.maximum(np.abs(xs[mask]), np.abs(ys[mask]))
    return np.mean(e)

def max_abs_err(xs, ys):
    e = np.abs(xs - ys)
    e = e[e < 90000]
    return np.max(e)


def mean_abs_err(xs, ys):
    e = np.abs(xs - ys)
    # Remove wet/dry outliers
    e = e[e < 90000]
    return np.mean(e)

def last_snapshot(filename):
    f = pd.read_csv(filename,
                    sep=r"\s+",
                    lineterminator='\n',
                    skiprows=0)
    num_nodes = int(f.iloc[0,1])
    last = f.iloc[-num_nodes:, 1].to_numpy()
    return last, num_nodes

def run_serial(path, rtol=0.05, atol=0.01):
    try:
        subprocess.run(os.path.join(WORK, "dgswem_serial"), check=True, cwd=path)
    except FileNotFoundError:
        pytest.skip("dgswem_serial executable not found...")

    d1, _ = last_snapshot(os.path.join(path , "fort.63.true"))
    d2, _ = last_snapshot(os.path.join(path , "fort.63"))
    np.testing.assert_allclose(d2, d1, rtol=rtol, atol=atol)

def run_parallel(path, rtol=0.05, atol=0.01):
    if not os.path.exists(os.path.join(WORK, "dgswem")):
        pytest.skip("dgswem executable not found...")

    subprocess.run(os.path.join(WORK, "adcprep") + " < in.prep", check=True, cwd=path, shell=True)
    subprocess.run("mpirun -np 2 " + os.path.join(WORK, "dgswem"), check=True, cwd=path, shell=True)
    subprocess.run(os.path.join(WORK, "adcpost") + " < out.prep", check=True, cwd=path, shell=True)

    d1, _ = last_snapshot(os.path.join(path , "fort.63.true"))
    d2, _ = last_snapshot(os.path.join(path , "fort.63"))
    np.testing.assert_allclose(d2, d1, rtol=rtol, atol=atol)

