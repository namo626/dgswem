import numpy as np
import os
import pytest
from util import *

case = "wetdry"
ADCIRC = "../../cpu/dgswem_serial"
ADCIRC_CUDA = "../../work/dgswem_serial"


@pytest.fixture(scope="module")
def adg_solution():
    return make_solution(case, ADCIRC, "serial")

@pytest.fixture(scope="module")
def cuda_solution():
    return make_solution(case, ADCIRC_CUDA, "cuda")


def test_cuda(adg_solution, cuda_solution):
    np.testing.assert_allclose(cuda_solution, adg_solution, atol=1e-3)
