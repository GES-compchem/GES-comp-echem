import pytest, pathlib
import numpy as np

from os import listdir
from os.path import abspath, dirname, join
from numpy.testing import assert_array_almost_equal, assert_almost_equal

from compechem.core.geometry import MolecularGeometry

# Get the path of the tests directory
TEST_DIR = dirname(abspath(__file__))

# Test the MolecularGeometry class constructor under normal conditions
def test_MolecularGeometry___init__():
    try:
        MolecularGeometry()
    except:
        assert False, "Exception raised during MolecularGeometry object construction"
    else:
        assert True

# Test the MolecularGeometry class from_xyz classmethod
def test_MolecularGeometry_from_xyz():
    
    xyzfile = join(TEST_DIR, "utils/xyz_examples/with_comment.xyz")

    try:
        _ = MolecularGeometry.from_xyz(xyzfile)
    except:
        assert False, "Exception raised during MolecularGeometry object construction"
    else:
        assert True

# Test the MolecularGeometry load_xyz method
def test_MolecularGeometry_load_xyz():

    folder = join(TEST_DIR, "utils/xyz_examples")
    for xyzfile in listdir(folder):
        print(xyzfile)

        path = join(folder, xyzfile)

        try:
            mol = MolecularGeometry()
            mol.load_xyz(path)
        
        except:
            assert False, f"Exception raised when loading the {xyzfile} file"

        expected = [
            np.array([-3.21035, -0.58504, -0.01395]),
            np.array([-2.24247, -0.61827, 0.01848]),
            np.array([-3.48920, -1.24911, 0.63429])
        ]

        assert mol.atomcount == 3
        assert mol.atoms == ["O", "H", "H"]
        assert_array_almost_equal(expected, mol.coordinates, decimal=6)

# Test the MolecularGeometry class __getitem__, __iter__ and __len__ methods
def test_MolecularGeometry_special_methods():
    
    xyzfile = join(TEST_DIR, "utils/xyz_examples/with_comment.xyz")
    mol = MolecularGeometry.from_xyz(xyzfile)

    # Test the len method
    assert len(mol) == 3

    # Test the getitem method
    atom, coord = mol[1]
    assert atom == "H"
    assert_array_almost_equal(coord, np.array([-2.24247, -0.61827, 0.01848]), decimal=6)

    # Test the failure of the getitem method when calling an invalid index
    try:
        _, _ = mol[3]
    except:
        assert True
    else:
        assert False, "An exception was expected when trying to access index 3"
    
    try:
        _, _ = mol[-1]
    except:
        assert True
    else:
        assert False, "An exception was expected when trying to access index -1"
    
    # Test the iterator method
    expected_atoms = ["O", "H", "H"]
    expected_coordinates = [
        np.array([-3.21035, -0.58504, -0.01395]),
        np.array([-2.24247, -0.61827, 0.01848]),
        np.array([-3.48920, -1.24911, 0.63429])
    ]

    for i, (atom, coord) in enumerate(mol):
        assert expected_atoms[i] == atom
        assert_array_almost_equal(expected_coordinates[i], coord, decimal=6)

# Test the append method
def test_MolecularGeometry_append():

    xyzfile = join(TEST_DIR, "utils/xyz_examples/with_comment.xyz")
    mol = MolecularGeometry.from_xyz(xyzfile)

    mol.append("Am", [0., 1., 2.])

    assert len(mol) == 4
    assert mol.atomcount == 4
    assert mol.atoms == ["O", "H", "H", "Am"]
    assert mol[3][0] == "Am"
    assert_array_almost_equal(mol[3][1], [0, 1, 2], decimal=6)


# Test the write_xyz method
def test_MolecularGeometry_write_xyz(tmp_path_factory):
    
    xyzfile = join(TEST_DIR, "utils/xyz_examples/with_comment.xyz")
    mol = MolecularGeometry.from_xyz(xyzfile)
    mol.append("N", np.array([0, 0, 0]))

    tmpdir = tmp_path_factory.mktemp("temp_xyz")
    outfile = join(tmpdir, "newxyz.xyz")

    mol.write_xyz(outfile, comment="This is a new xyz file")

    with open(outfile, "r") as file:
        lines = file.readlines()

    assert len(lines) == 6
    assert lines[0] == "4\n"
    assert lines[1] == "This is a new xyz file\n"
    assert lines[2] == "O\t-3.2103500000\t-0.5850400000\t-0.0139500000\n"
    assert lines[3] == "H\t-2.2424700000\t-0.6182700000\t0.0184800000\n"
    assert lines[4] == "H\t-3.4892000000\t-1.2491100000\t0.6342900000\n"
    assert lines[5] == "N\t0.0000000000\t0.0000000000\t0.0000000000\n"


# Test the remaining MolecularGeometry class properties
def test_MolecularGeometry_properties():

    xyzfile = join(TEST_DIR, "utils/xyz_examples/with_comment.xyz")
    mol = MolecularGeometry.from_xyz(xyzfile)

    assert_almost_equal(mol.mass, 18.01528, decimal=4)