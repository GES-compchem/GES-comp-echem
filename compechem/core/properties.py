from typing import Dict, List, Tuple, Union

from compechem.core.basewrapper import BaseWrapper


class Properties:
    """
    Class containing the properties associated to a system using a given level of theory.
    
    Attributes
    ----------
    electronic_energy: float
        Electronic energy in Hartree.
    vibronic_energy: float
        Vibronic energy in Hartree.
    pka: float
        The computed value of the pKa.
    mulliken_charges: List[float]
        The list of Mulliken charges.
    mulliken_spin_populations: List[float]
        The list of Mulliken spin populations.
    condensed_fukui_mulliken: List[float]
        The list of Fukui condensed functions computed from the Mulliken charges.
    """

    def __init__(self):
        self.electronic_energy: float = None
        self.vibronic_energy: float = None
        self.pka: float = None
        self.mulliken_charges: List[float] = []
        self.mulliken_spin_populations: List[float] = []
        self.condensed_fukui_mulliken: List[float] = []
    
    def __str__(self) -> str:
        ## IMPLEMENT IN MORE DETAILS
        msg = "--- Energies (Eh) ---\n"
        msg += f"Electronic energy: {self.electronic_energy} Eh\n"
        msg += f"Vibronic energy: {self.vibronic_energy} Eh\n"
        return msg



class PropertiesArchive:
    """
    Simple wrapper around a dictionary capable of automatically storing the properties of a
    given system. Once a wrapper is provided the level of theory is stored and the computed
    properties are automatically associated to the correct entry. 
    """
    def __init__(self) -> None:
        self.__properties: Dict[str, Properties] = {}
    
    def __get_key(self, wrapper: Union[str, BaseWrapper]) -> str:
        """
        Generates a univocal key capable of identifying the level of theory at which a set of
        properties have been computed.

        Arguments
        ---------
        wrapper: Union[str, BaseWrapper]
            A string or a wrapper object indicating the level of theory to be used in accessing
            the dictionary entries.
        
        Raises
        ------
        TypeError
            Exception raised if the wrapper argument does not match the expected types.
        """
        if type(wrapper) != str and not isinstance(wrapper, BaseWrapper):
            raise TypeError("The wrapper type must be str or be derived from BaseWrapper")

        return wrapper if type(wrapper) == str else wrapper.wrapper_id

    def __getitem__(self, wrapper: Union[str, BaseWrapper]) -> Properties:
        """
        Access the properties stored under a given level of theory.

        Arguments
        ---------
        wrapper: Union[str, BaseWrapper]
            A string or a wrapper object indicating the level of theory to be used in accessing
            the dictionary entries.
        
        Raises
        ------
        ValueError
            Exceprion raised if the specified wrapper id is not available in the archive.
        
        Returns
        -------
        Properties
            The Properties object associated with the selected level of theory.
        """
        key = self.__get_key(wrapper)

        if key not in self.__properties:
            raise ValueError(f"Cannot find any property for the level of theory {key}.")
        else:
            return self.__properties[key]
    
    def __setitem__(self, wrapper: Union[str, BaseWrapper], properties: Properties) -> None:
        """
        Sets the properties associated to a given level of theory. If the key does not exist 
        a new entry will be created.

        Arguments
        ---------
        wrapper: Union[str, BaseWrapper]
            A string or a wrapper object indicating the level of theory used in generating the
            properties.
        properties: Properties
            The properties to be associated to the specified level of theory.
        
        Raises
        ------
        TypeError
            Exceprion raised if the given properties object is not of type `Properites`.
        """
        if type(properties) != Properties:
            raise TypeError("The properties argument must be of type Properties")

        key = self.__get_key(wrapper)
        self.__properties[key] = properties
    
    def __iter__(self) -> Tuple[str, Properties]:
        """
        Class iterator returning the sequece of levels of theory and associated properties.

        Yields
        ------
        str
            The wrapper id representing the level of theory associated to the properties.
        Properties
            The properties computed at the specified level of theory.
        """
        for wrapper_id, properties in self.__properties.items():
            yield wrapper_id, properties
    
    def add(self, wrapper: Union[str, BaseWrapper]) -> None:
        """
        When called adds, if not already existing, a level of theory in the properties dictionary.

        wrapper: Union[str, BaseWrapper]
            A string or a wrapper object indicating the level of theory used in generating the
            properties.
        """
        key = self.__get_key(wrapper)
        if key not in self.__properties:
            self.__properties[key] = Properties()
    
    def remove(self, wrapper: Union[str, BaseWrapper]) -> None:
        """
        Removes the specified level of theory from the properties dictionary.

        wrapper: Union[str, BaseWrapper]
            A string or a wrapper object indicating the level of theory used in generating the
            properties.
        """
        key = self.__get_key(wrapper)
        del self.__properties[key]