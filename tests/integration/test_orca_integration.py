import pytest

from compechem.engines.orca import OrcaInput
from compechem.systems import System, Ensemble

from os import listdir
from os.path import dirname, abspath, isfile
from shutil import rmtree

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_almost_equal

# Get the path of the tests directory
TEST_DIR = dirname(abspath(__file__))

# Test the OrcaInput class constructor
def test_OrcaInput___init__():

    try:
        engine = OrcaInput(
            method="HF", basis_set="def2-SVP", aux_basis="def2/J", solvent="water"
        )

    except:
        assert False, "Unenxpected exception raised during OrcaInput class construction"

    else:
        assert engine.method == "HF"
        assert (
            engine.level_of_theory
            == "OrcaInput || method: HF | basis: def2-SVP | solvent: water"
        )


# Test the spe() function on a radical cation water molecule in DMSO
def test_OrcaInput_spe():

    engine = OrcaInput(
        method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent="DMSO"
    )
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz", charge=1, spin=2)

    try:
        engine.spe(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during SPE calculation"

    else:
        assert mol.properties.level_of_theory_electronic == engine.level_of_theory
        assert_array_almost_equal(
            mol.properties.electronic_energy, -75.942595825106, decimal=6
        )

        expected_mulliken_charges = np.array([0.377902, 0.311149, 0.310950])
        assert_array_almost_equal(
            expected_mulliken_charges, mol.properties.mulliken_charges, decimal=4
        )
        expected_mulliken_spin_populations = np.array([1.044851, -0.022422, -0.022429])
        assert_array_almost_equal(
            expected_mulliken_spin_populations,
            mol.properties.mulliken_spin_populations,
            decimal=4,
        )

        rmtree("output_files")


# Test the spe() function on a radical cation water molecule in DMSO without the inplace option
def test_OrcaInput_spe_no_inplace():

    engine = OrcaInput(
        method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent="DMSO"
    )
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz", charge=1, spin=2)

    try:
        newmol = engine.spe(mol, ncores=4)
    except:
        assert False, "Unexpected exception raised during SPE calculation"

    else:
        assert newmol.properties.level_of_theory_electronic == engine.level_of_theory

        assert_array_almost_equal(
            newmol.properties.electronic_energy, -75.942595825106, decimal=6
        )

        expected_mulliken_charges = np.array([0.377902, 0.311149, 0.310950])
        assert_array_almost_equal(
            expected_mulliken_charges,
            newmol.properties.mulliken_charges,
            decimal=4,
        )

        expected_mulliken_spin_populations = np.array([1.044851, -0.022422, -0.022429])
        assert_array_almost_equal(
            expected_mulliken_spin_populations,
            newmol.properties.mulliken_spin_populations,
            decimal=4,
        )

        rmtree("output_files")


# Test the spe() function on a radical cation water molecule in DMSO without the inplace option
def test_OrcaInput_spe_CCSD():

    engine = OrcaInput(method="DLPNO-CCSD", basis_set="def2-SVP", aux_basis="AutoAux")
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz", charge=1, spin=2)

    try:
        newmol = engine.spe(mol, ncores=4)
    except:
        assert False, "Unexpected exception raised during SPE calculation"

    else:
        assert newmol.properties.level_of_theory_electronic == engine.level_of_theory

        assert_array_almost_equal(
            newmol.properties.electronic_energy, -75.731114338261, decimal=6
        )

        expected_mulliken_charges = np.array([0.391458, 0.304269, 0.304274])
        assert_array_almost_equal(
            expected_mulliken_charges,
            newmol.properties.mulliken_charges,
            decimal=4,
        )

        expected_mulliken_spin_populations = np.array([1.061085, -0.030540, -0.030544])
        assert_array_almost_equal(
            expected_mulliken_spin_populations,
            newmol.properties.mulliken_spin_populations,
            decimal=4,
        )

        rmtree("output_files")


# Test that the correct suffix is generated when forbidden symbol is used
def test_OrcaInput_forbidden():

    engine = OrcaInput(method="DLPNO-CCSD(T)", basis_set="6-311++G**", aux_basis="AutoAux")
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz", charge=1, spin=2)

    try:
        engine.spe(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised"

    assert (
        isfile("./output_files/water_1_2_orca_DLPNO-CCSD-T-_6-311++G--_vacuum_spe.out")
        == True
    ), "Output file not found"


# Test the opt() function on a water molecule in vacuum
def test_OrcaInput_opt():

    engine = OrcaInput(method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        engine.opt(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during geometry optimization"

    else:
        assert mol.properties.level_of_theory_electronic == engine.level_of_theory
        assert mol.properties.level_of_theory_vibronic == engine.level_of_theory

        assert_almost_equal(mol.properties.electronic_energy, -76.272686998306, decimal=6)
        assert_almost_equal(mol.properties.vibronic_energy, 0.00301009, decimal=6)
        assert_almost_equal(mol.properties.gibbs_free_energy, -76.26967691, decimal=6)

        expected_mulliken_charges = np.array([-0.285546, 0.142772, 0.142773])
        assert_array_almost_equal(
            expected_mulliken_charges, mol.properties.mulliken_charges, decimal=4
        )

        expected_geometry = [
            np.array([-3.216653, -0.578663, -0.020175]),
            np.array([-2.244047, -0.623851, 0.023928]),
            np.array([-3.481320, -1.249905, 0.635066]),
        ]
        assert_array_almost_equal(expected_geometry, mol.geometry.coordinates, decimal=6)

        rmtree("output_files")


# Test the freq() function on a water molecule in vacuum
def test_OrcaInput_freq():

    engine = OrcaInput(method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        engine.freq(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during frequency analysis"

    else:
        assert mol.properties.level_of_theory_electronic == engine.level_of_theory
        assert mol.properties.level_of_theory_vibronic == engine.level_of_theory

        assert_almost_equal(mol.properties.electronic_energy, -76.272562753586, decimal=6)
        assert_almost_equal(mol.properties.vibronic_energy, 0.00327855, decimal=6)
        assert_almost_equal(mol.properties.gibbs_free_energy, -76.26928420, decimal=6)

        expected_mulliken_charges = np.array([-0.285593, 0.142795, 0.142798])
        assert_array_almost_equal(
            expected_mulliken_charges, mol.properties.mulliken_charges, decimal=4
        )

        rmtree("output_files")


# Test the nfreq() function on a water molecule in ethanol
def test_OrcaInput_nfreq():

    engine = OrcaInput(
        method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent="ethanol"
    )
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        engine.freq(mol, ncores=4, inplace=True)
    except:
        assert False, "Unexpected exception raised during numerical frequency analysis"

    else:
        assert mol.properties.level_of_theory_electronic == engine.level_of_theory
        assert mol.properties.level_of_theory_vibronic == engine.level_of_theory

        assert_almost_equal(mol.properties.electronic_energy, -76.283368184519, decimal=6)
        assert_almost_equal(mol.properties.vibronic_energy, 0.00317889, decimal=6)
        assert_almost_equal(mol.properties.gibbs_free_energy, -76.28018929, decimal=6)

        expected_mulliken_charges = np.array([-0.363774, 0.181883, 0.181890])
        assert_array_almost_equal(
            expected_mulliken_charges, mol.properties.mulliken_charges, decimal=4
        )

        rmtree("output_files")


# Test the scan() function on a water molecule in vacuum
def test_OrcaInput_scan():

    engine = OrcaInput(method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        ensemble: Ensemble = engine.scan(mol, scan="B 0 1 = 0.8, 1.5, 10", ncores=4)
    except:
        assert False, "Unexpected exception raised during relaxed surface scan"

    else:
        assert len(ensemble.systems) == 10

        expected_energies = np.array(
            [
                -76.23067389,
                -76.26216326,
                -76.27234942,
                -76.27002128,
                -76.26050938,
                -76.24706330,
                -76.23167352,
                -76.21556679,
                -76.19949798,
                -76.18392213,
            ]
        )

        calculated_energies = np.array(
            [system.properties.electronic_energy for system in ensemble.systems]
        )

        assert_array_almost_equal(calculated_energies, expected_energies, decimal=6)

        rmtree("output_files")


# Test the catching of runtime errors
def test_OrcaInput_runtime_error_input():

    engine = OrcaInput(method="PBU", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        engine.spe(mol, ncores=4)
    except:
        assert True
    else:
        assert False, "An exception was not raised on wrong input file."

    for filename in listdir("./"):
        if filename.endswith("_spe"):
            rmtree(filename)


def test_OrcaInput_runtime_error_missing_basis():

    engine = OrcaInput(
        method="DLPNO-CCSD", basis_set="def2-SVP", aux_basis="def2/J", solvent=None
    )
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz")

    try:
        engine.spe(mol, ncores=4)
    except:
        assert True
    else:
        assert False, "An exception was not raised on missing basis-set."

    for filename in listdir("./"):
        if filename.endswith("_spe"):
            rmtree(filename)


def test_OrcaInput_runtime_error_scf_not_converged():

    engine = OrcaInput(
        method="PBE",
        basis_set="def2-SVP",
        aux_basis="def2/J",
        solvent=None,
        scf_block={"maxiter": 2},
    )
    mol = System(f"{TEST_DIR}/utils/xyz_files/europium-aquoion.xyz", charge=3, spin=7)

    try:
        engine.spe(mol, ncores=4)
    except:
        assert True
    else:
        assert False, "An exception was not raised on SCF not converged."

    for filename in listdir("./"):
        if filename.endswith("_spe"):
            rmtree(filename)


def test_OrcaInput_runtime_error_wrong_multiplicity():

    engine = OrcaInput(method="PBE", basis_set="def2-SVP", aux_basis="def2/J", solvent=None)
    mol = System(f"{TEST_DIR}/utils/xyz_files/water.xyz", charge=0, spin=2)

    try:
        engine.spe(mol, ncores=4)
    except:
        assert True
    else:
        assert False, "An exception was not raised on wrong multiplicity."

    for filename in listdir("./"):
        if filename.endswith("_spe"):
            rmtree(filename)
