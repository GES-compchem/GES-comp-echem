import numpy as np

from typing import List, Tuple


class VibrationalData:

    def __init__(self) -> None:
        self.frequencies: List[float] = []
        self.normal_modes: List[np.ndarray] = []
        self.ir_transitions: List[Tuple[int, float]] = []
        self.ir_combination_bands: List[Tuple[int, int, float]] = []
        self.raman_transitions: List[Tuple[int, float, float]] = []