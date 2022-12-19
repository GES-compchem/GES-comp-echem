import logging, warnings
import compechem.config

from typing import Dict, List
from compechem.core.base import BaseEngine



logger = logging.getLogger(__name__)


class Properties:
    """
    Class containing the properties associated to a system using a given level of theory.

    Attributes
    ----------
    electronic_energy: float
        Electronic energy in Hartree.
    vibronic_energy: float
        Vibronic contribution to the total energy in Hartree.
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
        self.__level_of_theory_electronic: str = None
        self.__level_of_theory_vibronic: str = None

        self.__electronic_energy: float = None
        self.__vibronic_energy: float = None
        self.__helmholtz_free_energy: float = None
        self.__gibbs_free_energy: float = None
        self.__pka: float = None
        self.__mulliken_charges: List[float] = []
        self.__mulliken_spin_populations: List[float] = []
        self.__condensed_fukui_mulliken: Dict[str, List[float]] = {}
        self.__hirshfeld_charges: List[float] = []
        self.__hirshfeld_spin_populations: List[float] = []
        self.__condensed_fukui_hirshfeld: Dict[str, List[float]] = {}

    def __clear_electronic(self):
        self.__level_of_theory_electronic = None
        self.__electronic_energy = None
        self.__helmholtz_free_energy = None
        self.__gibbs_free_energy = None
        self.__pka = None
        self.__mulliken_charges = []
        self.__mulliken_spin_populations = []
        self.__condensed_fukui_mulliken = {}
        self.__hirshfeld_charges = []
        self.__hirshfeld_spin_populations = []
        self.__condensed_fukui_hirshfeld = {}

    def __clear_vibronic(self):
        self.__level_of_theory_vibronic = None
        self.__vibronic_energy = None
        self.__helmholtz_free_energy = None
        self.__gibbs_free_energy = None
        self.__pka = None
    
    def __check_engine(self, engine: BaseEngine) -> None:
        if not isinstance(engine, BaseEngine):
            raise TypeError("The engine argument must be derived from `BaseEngine`")

    def __validate_electronic(self, engine: BaseEngine) -> None:

        self.__check_engine(engine)

        if self.__level_of_theory_electronic is None:
            self.__level_of_theory_electronic = engine.level_of_theory

        elif self.__level_of_theory_electronic != engine.level_of_theory:
            if compechem.config.STRICT_MODE == True:
                msg = "Different electronic levels of theory used for calculating properties. Clearing properties with different electronic level of theory."
                logger.warning(msg)
                self.__clear_electronic()
                self.__level_of_theory_electronic = engine.level_of_theory

            else:
                msg = "Different electronic levels of theory used for calculating properties. Setting level of theory to undefined."
                logger.warning(msg)
                warnings.warn(msg)
                self.__level_of_theory_electronic = "Undefined"

    def __validate_vibronic(self, engine: BaseEngine) -> None:

        self.__check_engine(engine)

        if self.__level_of_theory_vibronic is None:
            self.__level_of_theory_vibronic = engine.level_of_theory
            
            if self.__pka is not None:
                if compechem.config.STRICT_MODE == True:
                    msg = "Added vibronic energy. Clearing pKa computed with electronic energy only."
                    logger.warning(msg)
                    self.__pka = None
                else:
                    msg = "Added vibronic energy to Properties with pKa previously computed with electronic energy only."
                    logger.warning(msg)
                    warnings.warn(msg)

        elif self.__level_of_theory_vibronic != engine.level_of_theory:

            if compechem.config.STRICT_MODE == True:
                msg = "Different vibronic levels of theory used for calculating properties. Clearing properties with different vibronic level of theory."
                logger.warning(msg)
                self.__clear_vibronic()
                self.__level_of_theory_vibronic = engine.level_of_theory

            else:
                msg = "Different vibronic levels of theory used for calculating properties. Setting level of theory to undefined."
                logger.warning(msg)
                warnings.warn(msg)
                self.__level_of_theory_vibronic = "Undefined"

    @property
    def level_of_theory_electronic(self) -> str:
        return self.__level_of_theory_electronic

    @property
    def level_of_theory_vibronic(self) -> str:
        return self.__level_of_theory_vibronic

    @property
    def electronic_energy(self) -> float:
        return self.__electronic_energy

    def set_electronic_energy(self, value: float, electronic_engine: BaseEngine) -> None:
        self.__validate_electronic(electronic_engine)
        self.__electronic_energy = value

    @property
    def vibronic_energy(self) -> float:
        return self.__vibronic_energy

    def set_vibronic_energy(self, value: float, vibronic_engine: BaseEngine) -> None:
        self.__validate_vibronic(vibronic_engine)
        self.__vibronic_energy = value

    @property
    def helmholtz_free_energy(self) -> float:
        return self.__helmholtz_free_energy

    def set_helmholtz_free_energy(
        self, value: float, electronic_engine: BaseEngine, vibronic_engine: BaseEngine
    ) -> float:
        self.__validate_electronic(electronic_engine)
        self.__validate_vibronic(vibronic_engine)
        self.__helmholtz_free_energy = value

    @property
    def gibbs_free_energy(self) -> float:
        return self.__gibbs_free_energy

    def set_gibbs_free_energy(
        self, value: float, electronic_engine: BaseEngine, vibronic_engine: BaseEngine
    ) -> float:
        self.__validate_electronic(electronic_engine)
        self.__validate_vibronic(vibronic_engine)
        self.__gibbs_free_energy = value

    @property
    def pka(self) -> float:
        return self.__pka

    def set_pka(
        self,
        value: float,
        electronic_engine: BaseEngine,
        vibronic_engine: BaseEngine = None,
    ) -> float:
        self.__validate_electronic(electronic_engine)
        if vibronic_engine is not None:
            self.__validate_vibronic(vibronic_engine)
        self.__pka = value

    @property
    def mulliken_charges(self) -> List[float]:
        return self.__mulliken_charges

    def set_mulliken_charges(
        self, value: List[float], electronic_engine: BaseEngine
    ) -> None:
        self.__validate_electronic(electronic_engine)
        self.__mulliken_charges = value

    @property
    def mulliken_spin_populations(self) -> List[float]:
        return self.__mulliken_spin_populations

    def set_mulliken_spin_populations(
        self, value: List[float], electronic_engine: BaseEngine
    ) -> None:
        self.__validate_electronic(electronic_engine)
        self.__mulliken_spin_populations = value

    @property
    def condensed_fukui_mulliken(self) -> Dict[str, List[float]]:
        return self.__condensed_fukui_mulliken

    def set_condensed_fukui_mulliken(
        self, value: Dict[str, List[float]], electronic_engine: BaseEngine
    ) -> None:
        self.__validate_electronic(electronic_engine)
        self.__condensed_fukui_mulliken = value
    
    @property
    def hirshfeld_charges(self) -> List[float]:
        return self.__hirshfeld_charges
    
    def set_hirshfeld_charges(
        self, value: List[float], electronic_engine: BaseEngine
    ) -> None:
        self.__validate_electronic(electronic_engine)
        self.__hirshfeld_charges = value
    
    @property
    def hirshfeld_spin_populations(self) -> List[float]:
        return self.__hirshfeld_spin_populations
    
    def set_hirshfeld_spin_populations(
        self, value: List[float], electronic_engine: BaseEngine
    ) -> None:
        self.__validate_electronic(electronic_engine)
        self.__hirshfeld_spin_populations = value
    
    @property
    def condensed_fukui_hirshfeld(self) -> Dict[str, List[float]]:
        return self.__condensed_fukui_hirshfeld

    def set_condensed_fukui_hirshfeld(
        self, value: Dict[str, List[float]], electronic_engine: BaseEngine
    ) -> None:
        self.__validate_electronic(electronic_engine)
        self.__condensed_fukui_hirshfeld = value
    

