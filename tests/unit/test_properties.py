import pytest
from numpy.testing import assert_array_almost_equal

import compechem.config as cc
from compechem.core.base import BaseEngine
from compechem.core.properties import Properties


def test_Properties___init__():

    try:
        Properties()

    except:
        assert False, "Unexpected exception raised on class construction"

    else:
        assert True


def test_Properties_properties():

    p = Properties()

    # Check that all the properties are empty on construction
    assert p.level_of_theory_electronic == None
    assert p.level_of_theory_vibronic == None
    assert p.electronic_energy == None
    assert p.vibronic_energy == None
    assert p.helmholtz_free_energy == None
    assert p.gibbs_free_energy == None
    assert p.pka == None
    assert p.mulliken_charges == {}
    assert p.mulliken_spin_populations == {}
    assert p.condensed_fukui_mulliken == {}

    # Define a BaseEngine instance to be used in setting the level of theory
    el_engine = BaseEngine("ElMethod")
    vib_engine = BaseEngine("VibMethod")

    # Set all properties
    p.set_electronic_energy(1, el_engine)
    p.set_vibronic_energy(2, vib_engine)
    p.set_helmholtz_free_energy(3, el_engine, vib_engine)
    p.set_gibbs_free_energy(4, el_engine, vib_engine)
    p.set_pka(5, el_engine, vib_engine)
    p.set_mulliken_charges({1: 6, 2: 7, 3: 8}, el_engine)
    p.set_mulliken_spin_populations({1: 9, 2: 10, 3: 11}, el_engine)
    p.set_condensed_fukui_mulliken({"f+": [0, 1, 2]}, el_engine)

    # Check that all the properties matces the set values
    assert p.level_of_theory_electronic == el_engine.level_of_theory
    assert p.level_of_theory_vibronic == vib_engine.level_of_theory
    assert p.electronic_energy == 1
    assert p.vibronic_energy == 2
    assert p.helmholtz_free_energy == 3
    assert p.gibbs_free_energy == 4
    assert p.pka == 5
    
    assert_array_almost_equal(
        [x for x in p.mulliken_charges.values()], [6, 7, 8], decimal=6
    )

    assert_array_almost_equal(
        [x for x in p.mulliken_spin_populations.values()], [9, 10, 11], decimal=6
    )

    assert_array_almost_equal(p.condensed_fukui_mulliken["f+"], [0, 1, 2], decimal=6)


def test_strict_mode_conflict_electronic():

    cc.STRICT_MODE = True

    p = Properties()
    first = BaseEngine("FirstMethod")
    second = BaseEngine("SecondMethod")

    assert first.level_of_theory != second.level_of_theory

    p.set_electronic_energy(0.1, first)
    p.set_mulliken_charges({1: 1, 2: 2, 3: 3}, second)

    assert p.electronic_energy == None
    assert p.level_of_theory_electronic == second.level_of_theory
    assert p.level_of_theory_vibronic == None
    assert_array_almost_equal(
        [x for x in p.mulliken_charges.values()], [1, 2, 3], decimal=6
    )


@pytest.mark.filterwarnings("ignore")
def test_not_strict_mode_conflict_electronic():

    cc.STRICT_MODE = False

    p = Properties()
    first = BaseEngine("FirstMethod")
    second = BaseEngine("SecondMethod")

    assert first.level_of_theory != second.level_of_theory

    p.set_electronic_energy(0.1, first)
    p.set_mulliken_charges({1: 1, 2: 2, 3: 3}, second)

    assert p.electronic_energy == 0.1
    assert p.level_of_theory_electronic == "Undefined"
    assert p.level_of_theory_vibronic == None
    assert_array_almost_equal(
        [x for x in p.mulliken_charges.values()], [1, 2, 3], decimal=6
    )


def test_strict_mode_conflict_vibronic():

    cc.STRICT_MODE = True

    p = Properties()
    first = BaseEngine("FirstMethod")
    second = BaseEngine("SecondMethod")

    assert first.level_of_theory != second.level_of_theory

    p.set_vibronic_energy(0.1, first)
    p.set_gibbs_free_energy(0.6, first, second)

    assert p.vibronic_energy == None
    assert p.level_of_theory_electronic == first.level_of_theory
    assert p.level_of_theory_vibronic == second.level_of_theory
    assert p.gibbs_free_energy == 0.6


@pytest.mark.filterwarnings("ignore")
def test_not_strict_mode_conflict_vibronic():

    cc.STRICT_MODE = False

    p = Properties()
    first = BaseEngine("FirstMethod")
    second = BaseEngine("SecondMethod")

    assert first.level_of_theory != second.level_of_theory

    p.set_vibronic_energy(0.1, first)
    p.set_gibbs_free_energy(0.6, first, second)

    assert p.vibronic_energy == 0.1
    assert p.level_of_theory_electronic == first.level_of_theory
    assert p.level_of_theory_vibronic == "Undefined"
    assert p.gibbs_free_energy == 0.6


def test_pka_vibronic_addition_strict():

    cc.STRICT_MODE = True

    p = Properties()
    first = BaseEngine("FirstMethod")
    second = BaseEngine("SecondMethod")

    assert first.level_of_theory != second.level_of_theory

    p.set_pka(0.0012, first)

    assert p.pka == 0.0012
    assert p.level_of_theory_electronic == first.level_of_theory
    assert p.level_of_theory_vibronic == None
    assert p.vibronic_energy == None

    p.set_vibronic_energy(1.0, second)

    assert p.pka == None
    assert p.level_of_theory_electronic == first.level_of_theory
    assert p.level_of_theory_vibronic == second.level_of_theory
    assert p.vibronic_energy == 1.0


@pytest.mark.filterwarnings("ignore")
def test_pka_vibronic_addition_not_strict():

    cc.STRICT_MODE = False

    p = Properties()
    first = BaseEngine("FirstMethod")
    second = BaseEngine("SecondMethod")

    assert first.level_of_theory != second.level_of_theory

    p.set_pka(0.0012, first)

    assert p.pka == 0.0012
    assert p.level_of_theory_electronic == first.level_of_theory
    assert p.level_of_theory_vibronic == None
    assert p.vibronic_energy == None

    p.set_vibronic_energy(1.0, second)

    assert p.pka == 0.0012
    assert p.level_of_theory_electronic == first.level_of_theory
    assert p.level_of_theory_vibronic == second.level_of_theory
    assert p.vibronic_energy == 1.0


def test_check_engine():

    p = Properties()

    try:
        p.set_electronic_energy(0.1, "This is a string")
    except:
        assert True
    else:
        assert False, "No exception raised when wrong type has been given as engine"

    try:
        p.set_vibronic_energy(0.5, 1)
    except:
        assert True
    else:
        assert False, "No exception raised when wrong type has been given as engine"
