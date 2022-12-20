import pytest

from compechem.engines.orca import OrcaInput
from compechem.systems import System
from os.path import dirname, abspath
from shutil import rmtree

import numpy as np
from numpy.testing import assert_array_almost_equal

# Get the path of the tests directory
TEST_DIR = dirname(abspath(__file__))

# Test the OrcaInput class constructor
def test_OrcaInput___init__():

    try:
        orca = OrcaInput(
            method="HF", basis_set="def2-SVP", aux_basis="def2/J", solvent="water"
        )

    except:
        assert False, "Unenxpected exception raised during OrcaInput class construction"

    else:
        assert orca.method == "HF"
        assert (
            orca.level_of_theory
            == "OrcaInput || method: HF | basis: def2-SVP | solvent: water"
        )


# Test the spe() function on a water molecule in vacuum
def test_OrcaInput_spe_vac():

    orca = OrcaInput(method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        orca.spe(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during SPE calculation"

    else:
        assert mol.properties.level_of_theory_electronic == orca.level_of_theory
        assert_array_almost_equal(
            mol.properties.electronic_energy, -76.272562753586, decimal=6
        )
        rmtree("output_files")


# Test the spe() function on a water molecule in solvent
def test_OrcaInput_spe_vac():

    orca = OrcaInput(
        method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent="water"
    )
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        orca.spe(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during SPE calculation"

    else:
        assert mol.properties.level_of_theory_electronic == orca.level_of_theory
        assert_array_almost_equal(
            mol.properties.electronic_energy, -76.283779209329, decimal=6
        )
        rmtree("output_files")


# Test the opt() function on a water molecule in vacuum
def test_OrcaInput_opt_vac():

    orca = OrcaInput(method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        orca.opt(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during geometry optimization"

    else:
        assert mol.properties.level_of_theory_electronic == orca.level_of_theory
        assert mol.properties.level_of_theory_vibronic == orca.level_of_theory

        assert_array_almost_equal(
            mol.properties.electronic_energy, -76.272686998306, decimal=6
        )
        assert_array_almost_equal(mol.properties.vibronic_energy, 0.00301009, decimal=6)

        expected = [
            np.array([-3.216653, -0.578663, -0.020175]),
            np.array([-2.244047, -0.623851, 0.023928]),
            np.array([-3.481320, -1.249905, 0.635066]),
        ]
        assert_array_almost_equal(expected, mol.geometry.coordinates, decimal=6)

        rmtree("output_files")
