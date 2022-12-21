import pytest

from compechem.engines.orca import OrcaInput

ORCADIR = "orcadir"

# Test the OrcaInput class constructor
def test_OrcaInput___init__():

    try:
        orca = OrcaInput(
            method="HF",
            basis_set="def2-SVP",
            aux_basis="def2/J",
            solvent="water",
            ORCADIR=ORCADIR,
        )

    except:
        assert False, "Unenxpected exception raised during OrcaInput class construction"

    else:
        assert orca.method == "HF"
        assert (
            orca.level_of_theory
            == "OrcaInput || method: HF | basis: def2-SVP | solvent: water"
        )
