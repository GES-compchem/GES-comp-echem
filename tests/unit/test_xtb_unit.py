import pytest

from compechem.engines.xtb import XtbInput

XTBPATH = "xtbpath"

# Test the XtbInput class constructor
def test_XtbInput___init__():

    try:
        xtb = XtbInput(
            method="gfn2",
            solvent="water",
            XTBPATH=XTBPATH,
        )

    except:
        assert False, "Unenxpected exception raised during XtbInput class construction"

    else:
        assert xtb.method == "gfn2"
        assert xtb.level_of_theory == "XtbInput || method: gfn2 | solvent: water"
