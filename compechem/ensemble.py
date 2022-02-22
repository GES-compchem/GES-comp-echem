import os
import numpy as np


class Ensemble:
    """Ensemble object, containing a series of Molecule objects.

    Attributes
    ----------
    name : str
        name of the molecule represented in the ensemble, taken from the first element of the 
        ensemble
    atomcount : int
        number of atoms contained each ensemble, taken from the first element of the ensemble
    energies : dict
        dictionary containing the electronic/vibronic energies of the molecule,
        calculated at various levels of theory
    """

    def __init__(self, molecules_list) -> None:
        """
        Parameters
        ----------
        molecules_list : list
            list of Molecule objects generating the ensemble.
        """

        self.name = molecules_list[0].name
        self.atomcount: molecules_list[0].atomcount

        self.molecules = [mol for mol in molecules_list]

        self.energies: dict = {}

    class Energies:
        """Ensemble energies.
        """

        def __init__(
            self, method: str = None, electronic: float = None, vibronic: float = None,
        ) -> None:
            """
            Parameters
            ----------
            method : str, optional
                level of theory, by default None
            electronic : float, optional
                electronic energy (in Hartree), by default None
            vibronic : float, optional
                vibronic contribution to the total energy (in Hartree),
                by default None
            """
            self.method = method
            self.electronic = electronic
            self.vibronic = vibronic

    def boltzmann_average(self, method_el, method_vib=None, temperature=297.15):
        """Calculates the average free Gibbs energy of the ensemble (in Hartree), weighted for each 
        molecule by its Boltzmann factor.
        
        Parameters
        ----------
        method_el : str
            level of theory for the electronic energies
        method_vib : str
            level of theory for the vibronic contributions, by default equal to method_el
        temperature : float
            temperature at which to calculate the Boltzmann average, by default 297.15 K
        """

        if method_vib is None:
            method_vib = method_el

        energies = []

        for mol in self.molecules:
            energies.append(mol.energies[method_el].electronic + mol.energies[method_vib].vibronic)

        # necessary for avoiding overflows with large energies
        relative_energies = [energy - min(energies) for energy in energies]

        boltzmann_constant = 3.167e-6  # Eh/K
        temperature = temperature

        partition_function = np.sum(
            np.exp([-energy / (boltzmann_constant * temperature) for energy in relative_energies])
        )

        populations = [
            np.exp(-energy / (boltzmann_constant * temperature)) / partition_function
            for energy in relative_energies
        ]

        weighted_energies = [
            energy * population for energy, population in zip(energies, populations)
        ]

        boltzmann_entropy = -boltzmann_constant * np.sum(populations * np.log(populations))

        self.energies[f"{method_el}"] = self.Energies(
            method=f"{method_el}",
            electronic=np.sum([weighted_energies]) - temperature * boltzmann_entropy,
            vibronic=0,
        )
