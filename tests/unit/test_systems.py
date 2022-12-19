import pytest
from os.path import abspath, dirname, join

from compechem.systems import System

# Get the path of the tests directory
TEST_DIR = dirname(abspath(__file__))

def test_System___init__():

    xyzfile = join(TEST_DIR, "utils/xyz_examples/with_comment.xyz")

    try:
        sys = System(xyzfile)
    
    except:
        assert False, "Exception raised during `System` class constructor"
    
    assert sys.geometry.atomcount == 3
    assert sys.charge == 0
    assert sys.spin == 1
    assert sys.box_side == None
    assert sys.is_periodic == False