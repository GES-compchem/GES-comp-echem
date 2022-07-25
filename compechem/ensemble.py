import os, shutil
import numpy as np
from compechem.molecule import Energies, System
import logging
from subprocess import getoutput
from itertools import chain

logger = logging.getLogger(__name__)


class Ensemble:
    """Ensemble object, containing a series of System objects.

    Attributes
    ----------
    name : str
        name of the system represented in the ensemble, taken from the first element of 
        the ensemble
    atomcount : int
        number of atoms contained each ensemble, taken from the first element of the 
        ensemble
    energies : dict
        dictionary containing the electronic/vibronic energies of the systems,
        calculated at various levels of theory
    """

    def __init__(self, systems_list) -> None:
        """
        Parameters
        ----------
        systems_list : any
            generator for creating the Ensemble.
            Options:
                - list of System objects
                - .xyz trajectory file
        """

        self.name = systems_list[0].name
        self.atomcount = systems_list[0].atomcount

        self.container = systems_list

        self.energies: dict = {}

    def __iter__(self):
        chained = chain.from_iterable(self.container)
        for item in chained:
            yield item

    def __getitem__(self, index):
        space = 0
        for iterator in self.container:
            space += len(iterator)
            if index < space:
                local_index = index - (space - len(iterator))
                return iterator[local_index]

    def __len__(self):
        return sum([len(iterator) for iterator in self.container])

    def add(self, iterator):
        """adds another file to the trajectory from an iterator object (e.g., MDTrajectory)

        Parameters
        ----------
        iterator : iterable
            iterator object (e.g., MDTrajectory object) 
        """
        self.container.append(iterator)

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
        """Calculates the average free Gibbs energy of the ensemble (in Hartree), weighted 
        for each molecule by its Boltzmann factor.

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

        NOTE: the vibronic contributions are included in the electronic component, which
        actually contains the TOTAL energy of the system. Maybe in the future I'll think of 
        how to separate the two contributions - LB
        """

        energies = []

        for sys in self.container:
            if method_vib is None:
                energies.append(sys.energies[method_el].electronic)
            else:
                energies.append(
                    sys.energies[method_el].electronic + sys.energies[method_vib].vibronic
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


class MDTrajectory:
    def __init__(self, system, method):
        """Iterator class for MD trajectories. Data is computed only when accessing the
        elements of the object (via __getitem__ or __iter__)

        Parameters
        ----------
        system : str
            name of the system outputted by an MD calculator
        method : str
            level of theory at which the simulation was ran
        """
        self.name = system
        self.method = method

        self.md_out = f"MD_data/{system}_md.out"
        self.geo_end = f"MD_data/{system}_geo_end.xyz"
        self.atomcount = int(getoutput(f"head -1 {self.geo_end}"))

        head = getoutput(f"head -50 {self.md_out}|grep step|head -2").split()
        mdrestartfreq = int(head[-1]) - int(head[-4])

        nframes = int(getoutput(f"tail -50 {self.md_out}|grep step|tail -1").split()[-1])

        self.frames = list(range(0, nframes + 1, mdrestartfreq))

    def __iter__(self):
        for index in self.frames:
            yield self[index]

    def __getitem__(self, index):
        """returns the System object corresponding to the requested index/MD step

        Parameters
        ----------
        index : int
            MD step in the simulation

        Returns
        -------
        system
            System object corresponding to the requested frame in the MD simulation
        """
        with open(self.geo_end, "r") as f:
            with open(f"{self.name}_{index}.xyz", "w+") as o:

                self.atomcount = int(f.readline())
                o.write(f"{self.atomcount}\n")
                start = None

                for i, line in enumerate(f):

                    if line == f"MD iter: {index}\n":
                        start = i
                        o.write(line)

                    if start and i > start and i < start + self.atomcount + 1:
                        o.write(
                            f"{line.split()[0]}\t{line.split()[1]}\t{line.split()[2]}\t{line.split()[3]}\n"
                        )

                    if start and i > start + self.atomcount + 1:
                        break

        system = System(f"{self.name}_{index}.xyz")
        os.remove(f"{self.name}_{index}.xyz")

        with open(self.md_out, "r") as f:
            found = False
            for line in f:
                if line == f"MD step: {index}\n":
                    found = True
                if "Total MD Energy" in line and found:
                    system.energies[self.method] = Energies(
                        method=self.method, electronic=line.split()[3], vibronic=None,
                    )
                    break

        return system

    def __len__(self):
        return len(self.frames)

