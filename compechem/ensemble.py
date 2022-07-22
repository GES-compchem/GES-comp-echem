import os
import numpy as np
from compechem.molecule import Energies, Molecule


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
        molecules_list : any
            generator for creating the Ensemble.
            Options:
                - list of Molecule objects
                - .xyz trajectory file
        """

        if type(molecules_list) == list and type(molecules_list[0] == Molecule):
            self.name = molecules_list[0].name
            self.atomcount = molecules_list[0].atomcount

            self.molecules = [mol for mol in molecules_list]

            self.energies: dict = {}

        elif type(molecules_list) == str and ".xyz" in molecules_list:

            self.name = os.path.basename(molecules_list).strip(".xyz")
            self.trajectory: str = molecules_list
            self.energies: dict = {}

            with open(self.trajectory, "r") as f:
                self.atomcount = int(f.readline())

    def read_energies(self, method):
        """reads energies from trajectory file (parsed by tools.save_dftb_trajectory()) and
        calculates the average energy for the given trajectory

        Parameters
        ----------
        method : str
            level of theory to assign to the self.energies dictionary key
        """

        energies = []

        with open(self.trajectory, "r") as f:
            for line in f:
                if "Energy" in line:
                    energies.append(float(line.split()[-2]))

        average_energy = sum(energies) / len(energies)

        self.energies[method] = Energies(
            method=method, electronic=average_energy, vibronic=0
        )

    def boltzmann_average(
        self, method_el: str, method_vib: str = None, temperature: float = 297.15
    ):
        """Calculates the average free Gibbs energy of the ensemble (in Hartree), weighted for each
        molecule by its Boltzmann factor.

        Parameters
        ----------
        method_el : str
            level of theory for the electronic energies
        method_vib : str, optional
            level of theory for the vibronic contributions, if method_vib is None (default),
            the ensemble energy will be evaluated only considering the electronic component
        temperature : float
            temperature at which to calculate the Boltzmann average, by default 297.15 K

        Returns
        -------
        Creates an Energies object with the total free Gibbs energy of the ensemble.

        NOTE: the vibronic contributions are included in the electronic component, which actually
        contains the TOTAL energy of the system. Maybe in the future I'll think of how to separate
        the two contributions - LB
        """

        energies = []

        for mol in self.molecules:
            if method_vib is None:
                energies.append(mol.energies[method_el].electronic)
            else:
                energies.append(
                    mol.energies[method_el].electronic + mol.energies[method_vib].vibronic
                )

        # necessary for avoiding overflows with large energies
        relative_energies = [energy - min(energies) for energy in energies]

        boltzmann_constant = 3.167e-6  # Eh/K
        temperature = temperature

        partition_function = np.sum(
            np.exp(
                [
                    -energy / (boltzmann_constant * temperature)
                    for energy in relative_energies
                ]
            )
        )

        populations = [
            np.exp(-energy / (boltzmann_constant * temperature)) / partition_function
            for energy in relative_energies
        ]

        weighted_energies = [
            energy * population for energy, population in zip(energies, populations)
        ]

        boltzmann_entropy = -boltzmann_constant * np.sum(populations * np.log(populations))

        self.energies[method_el] = Energies(
            method=method_el,
            electronic=np.sum([weighted_energies]) - temperature * boltzmann_entropy,
            vibronic=0,
        )
