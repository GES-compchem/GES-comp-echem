import pytest

from compechem.core.base import BaseEngine

# Test the BaseEngine class constructor
def test_BaseEngine___init__():

    try:
        obj = BaseEngine("DummyMethod")
    
    except:
        assert False, "Unenxpected exception raised during BaseEngine class construction"
    
    else:    
        assert obj.method == "DummyMethod"
        assert obj.level_of_theory == "BaseEngine || method: DummyMethod"


# Test that, when inheriting form BaseEngine, the expected behavior is obtained
def test_BaseEngine_inheritance():

    class DerivedEngine(BaseEngine):

        def __init__(self, method: str) -> None:
            super().__init__(method)
    
    try:
        dobj = DerivedEngine("DerivedMethod")
    
    except:
        assert False, "Unenxpected exception raised during inheritance from BaseEngine"
    
    else:    
        assert dobj.method == "DerivedMethod"
        assert dobj.level_of_theory == "DerivedEngine || method: DerivedMethod"