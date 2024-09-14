import numpy as np
import os
import pytest
from util import *

ADCIRC = "../../cpu/dgswem_serial"

ADCIRC_CUDA = "../../work/dgswem_serial"

@pytest.fixture(scope="module")
def adg_solution():
    try: 
        err = os.system("cd quarter_annular/ && %s --solver 1 > /dev/null && mv fort.63 fort.63.adg "% ADCIRC)
        if err != 0:
            raise Exception("DG-SWEM run error.")
    except:
        raise

    f14 = 'quarter_annular/fort.14'
    tr, d = read_63_all('quarter_annular/fort.63.adg', f14)
    total_nodes = len(tr.x)
    adg = d[-total_nodes:]
    return adg

@pytest.fixture(scope="module")
def cuda_solution():
    try:
        err = os.system("cd quarter_annular/ && %s --solver 1 > /dev/null && mv fort.63 fort.63.cuda "% ADCIRC_CUDA)
        if err != 0:
            raise Exception("DG-SWEM CUDA run error.")
    except:
        raise

    f14 = 'quarter_annular/fort.14'
    tr, d = read_63_all('quarter_annular/fort.63.cuda', f14)
    total_nodes = len(tr.x)
    adg = d[-total_nodes:]
    return adg

'''
@pytest.fixture(scope="module")
def mpi_solution():
    # MPI version
    os.system("cd quarter_annular/ && %s --np 2 --partmesh > /dev/null && %s --np 2 --prepall > /dev/null  " % (ADCPREP, ADCPREP))
    os.system("cd quarter_annular/ && mpirun -np 2 %s --solver 1 > /dev/null && mv fort.63 fort.63.par " % PADCIRC)

    f14 = 'quarter_annular/fort.14'
    tr, d = read_63_all('quarter_annular/fort.63.par', f14)
    total_nodes = len(tr.x)
    par = d[-total_nodes:]

    return par

@pytest.fixture(scope="module")
def mpi_cuda_solution():
    # MPI version
    os.system("cd quarter_annular/ && %s --np 2 --partmesh > /dev/null && %s --np 2 --prepall > /dev/null  " % (ADCPREP_CUDA, ADCPREP_CUDA))
    os.system("cd quarter_annular/ && mpirun -np 2 %s --solver 1 > /dev/null && mv fort.63 fort.63.par.cuda " % PADCIRC_CUDA)

    f14 = 'quarter_annular/fort.14'
    tr, d = read_63_all('quarter_annular/fort.63.par.cuda', f14)
    total_nodes = len(tr.x)
    par = d[-total_nodes:]

    return par

'''

def test_cuda(adg_solution, cuda_solution):
    np.testing.assert_allclose(cuda_solution, adg_solution, rtol=1e-8)

'''
def test_mpi_cuda(mpi_solution, mpi_cuda_solution):
    np.testing.assert_allclose(mpi_cuda_solution, mpi_solution, rtol=1e-8)
    '''
