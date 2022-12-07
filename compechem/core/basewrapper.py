from abc import ABC, abstractmethod

class BaseWrapper(ABC):
    """
    Simple abstract base class for the definition of a software wrapper. The class has a
    single abstract propery `wrapper_id` that, by default, automatically returns a string
    with the wrapper name and terminated by `: `.
    """
    def __init__(self) -> None:
        pass

    @property
    @abstractmethod
    def wrapper_id(self) -> str:
        return f"{self.__class__.__name__}: "
