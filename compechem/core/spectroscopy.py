import numpy as np
from typing import List, Tuple


class VibrationalData:
    """
    The `VibrationalData` class holds all the information related to the vibrational properties of a given system.

    Attributes
    ----------
    frequencies: List[float]
        The list of frequencies associated to each mode of vibration. The leading modes, associated to rotations and
        translations will be associated with zero frequencies.
    normal_modes: List[np.ndarray]
        The vectors encoding the cartesian displacements associated to each mode of vibration.
    ir_transitions: List[Tuple[int, float]]
        The list of tuples associated with each infrared transitions. Each tuple is composed by the index of the mode
        involved and the corresponding intensity in km/mol.
    ir_combination_bands: List[Tuple[int, int, float]]
        The list of tuples associated with each combination and overtone transitions. Each tuple is composed by the
        index of the two modes involved and the corresponding intensity in km/mol.
    raman_transitions: List[Tuple[int, float, float]]
        The list of tuples associated with each Raman transition. Each tuple is composed by the index of the mode
        involved, the corresponding activity and depolarization.
    """

    def __init__(self) -> None:
        self.frequencies: List[float] = []
        self.normal_modes: List[np.ndarray] = []
        self.ir_transitions: List[Tuple[int, float]] = []
        self.ir_combination_bands: List[Tuple[int, int, float]] = []
        self.raman_transitions: List[Tuple[int, float, float]] = []
    