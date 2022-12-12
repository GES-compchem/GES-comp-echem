from compechem.core.geometry import MolecularGeometry


# Test the MolecularGeometry class constructor under normal conditions
def test_MolecularGeometry___init__():
    try:
        MolecularGeometry()
    except:
        assert False, "Exception raised during MolecularGeometry object construction"
    else:
        assert True