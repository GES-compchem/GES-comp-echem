class BaseEngine:
    """
    Simple base class for the definition of a software wrapper. The class sets the 
    `wrapper_info` attribute initialized with the name of the wrapper itself terminated by `: `.

    Arguments
    ---------
    method: str
        The string indicating the method to be used in the calculation.
    
    """
    def __init__(self, method: str) -> None:
        self.method = method
        self.level_of_theory = f"{self.__class__.__name__} || method: {self.method}"
