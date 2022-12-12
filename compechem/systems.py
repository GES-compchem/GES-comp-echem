from __future__ import annotations

import os
import numpy as np
import logging

from typing import List

from compechem.constants import kB
from compechem.core.geometry import MolecularGeometry
from compechem.core.properties import PropertiesArchive

logger = logging.getLogger(__name__)


class System(MolecularGeometry):
    """
    The System object describes a generic molecular system, its geometry and its computed
    properties.

    Parameters
    ----------
    xyz_file : str
        path with the .xyz file containing the system geometry
    charge : int, optional
        total charge of the system. Defaults to 0 (neutral)
    spin : int, optional
        total spin of the system. Defaults to 1 (singlet)
    box_side : float, optional
        for periodic systems, defines the length (in Å) of the box side

    Attributes
    ----------
    name : str
        Name of the system, taken from the provided `.xyz` file
    charge : int
        Total charge of the system
    spin : int
        Total spin multeplicity of the system (2S+1)
    box_side : float
        For periodic systems, defines the length (in Å) of the box side
    properties : PropertiesArchive
        The PropertiesArchive object containing all the properties of the system computed
        using a specified level of theory.
    flags : list
        list containing all "warning" flags which might be encountered during calculations.
    """

    def __init__(
        self, xyz_file: str, charge: int = 0, spin: int = 1, box_side: float = None
    ) -> None:

        super().__init__()
        super().load_xyz(xyz_file)
        self.name = os.path.basename(xyz_file).strip(".xyz")
        self.charge: int = charge
        self.spin: int = spin
        self.box_side = box_side
        self.properties: PropertiesArchive = PropertiesArchive()
        self.flags: list = []

    def __str__(self):
        info = "=========================================================\n"
        info += f"SYSTEM: {self.name}\n"
        info += "=========================================================\n\n"
        info += f"Number of atoms: {self.atomcount}\n"
        info += f"Charge: {self.charge}\n"
        info += f"Spin multeplicity: {self.spin}\n"

        if self.box_side:
            info += f"Periodic system with box size: {self.box_side:.4f} Å\n"
        info += "\n"

        info += "********************** GEOMETRY *************************\n\n"
        info += f"Total system mass: {self.mass:.4f} amu\n\n"

        info += "----------------------------------------------\n"
        info += " index  atom    x (Å)      y (Å)      z (Å)   \n"
        info += "----------------------------------------------\n"
        for idx, (atom, coordinates) in enumerate(self):
            info += f" {idx:<6}{atom:^6}"
            for c in coordinates:
                info += "{0:^11}".format(f"{c:.5f}")
            info += "\n"
        info += "----------------------------------------------\n\n"

        if len(self.properties) != 0:
            info += "********************** PROPERTIES *************************\n\n"
            for level_of_theory, properties in self.properties:
                info += f"{level_of_theory}\n"
                info += str(properties)
                info += "\n"

        if self.flags != []:
            info += "********************** WARNINGS **************************\n\n"
            for warning in self.flags:
                info += f"{warning}\n"

        return info

    def load_xyz(self, xyz_file: str) -> None:
        """
        Updates the current geometry from an external `.xyz` file.

        Parameters
        ----------
        xyz_file : str
            path with the `.xyz` file of the geometry containing the new coordinates
        """
        super().load_xyz(xyz_file)
        self.properties.clear()

    def write_gen(self, gen_file: str, box_side: float = None):
        """Writes the current geometry to a .gen file.

        Parameters
        ----------
        gen_file : str
            path to the output `.gen` file
        box_side : float, optional
            for periodic systems, defines the length (in Å) of the box side
        """

        if box_side is None:
            box_side = self.box_side

        with open(gen_file, "w") as file:

            file.write(f" {str(self.atomcount)} ")
            file.write("S\n" if self.is_periodic else "C\n")

            atom_types = []
            for element in self.__atoms:
                if element not in atom_types:
                    atom_types.append(element)

            for atom in atom_types:
                file.write(f" {atom}")
            file.write("\n")

            i = 1
            for atom, coordinates in self:
                line = f"{atom}\t" + "\t".join(list(coordinates)) + "\n"
                for index, atom_type in enumerate(atom_types):
                    if line.split()[0] == atom_type:
                        file.write(f"{i} {line.replace(atom_type, str(index + 1))}")
                        i += 1

            if self.is_periodic:
                file.write(f" 0.000 0.000 0.000\n")
                file.write(f" {box_side} 0.000 0.000\n")
                file.write(f" 0.000 {box_side} 0.000\n")
                file.write(f" 0.000 0.000 {box_side}")

    @property
    def is_periodic(self) -> bool:
        """
        Indicates if the system is periodic or not

        Returns
        -------
        bool
            True if the system is periodic (the `box_side` is not None), False otherwise.
        """
        return True if self.box_side is not None else False


class Ensemble:
    """
    Ensemble object, containing a series of System objects.

    Parameters
    ----------
    systems : List[System]
        The list of System objects to be included in the Ensemble.

    Attributes
    ----------
    name : str
        Name of the system represented in the ensemble, taken from the first element of
        the ensemble
    systems : List[System]
        The list of System objects in the Ensemble.
    properties : PropertiesArchive
        The property archive containing the average ensamble properties calculated at
        various levels of theory.
    """
    def __init__(self, systems: List[System]) -> None:

        if len(systems) == 0:
            raise ValueError("Cannot operate on an empty systems array")

        if any(system != systems[0] for system in systems):
            raise RuntimeError("Different systems encountered in list")

        self.name: str = systems[0].name
        self.systems: List[System] = systems
        self.properties: PropertiesArchive = PropertiesArchive()

    def __iter__(self) -> System:
        for item in self.systems:
            yield item

    def __getitem__(self, index: int) -> System:
        if index < 0 or index >= len(self.systems):
            raise ValueError("Index out of bounds")

        return self.systems[index]

    def __len__(self) -> int:
        return len(self.systems)

    @property
    def atomcount(self) -> int:
        """
        The number of atoms in the system

        Returns
        -------
        int
            The total number of atoms
        """
        return self.systems[0].atomcount

    def add(self, systems: List[System]):
        """
        Append more Systems to the ensemble

        Parameters
        ----------
        systems : List[System]
            The list of systems to be added to the ensamble
        """
        if any(system != self.systems[0] for system in systems):
            raise RuntimeError("Different systems encountered in list")

        for system in systems:
            self.systems.append(system)

    def boltzmann_average(
        self,
        method_el: str,
        method_vib: str = None,
        temperature: float = 297.15,
    ):
        """
        Calculates the average free Gibbs energy of the ensemble (in Hartree), weighted
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

        for system in self.systems:
            if method_vib is None:
                energies.append(system.properties[method_el].electronic_energy)
            else:
                energies.append(
                    system.properties[method_el].electronic_energy
                    + system.properties[method_vib].vibronic_energy
                )

        # Compute the relative energy of each system in respect to the minimum to avoid overflows
        # when computing exponential of large magnitude values
        dE = [energy - min(energies) for energy in energies]

        # Compute the relative partition function starting from the relative energy list
        relative_Z = np.sum(np.exp([-energy / (kB * temperature) for energy in dE]))

        # Compute the populations for each system given the boltzmann distribution
        populations = [np.exp(-energy / (kB * temperature)) / relative_Z for energy in dE]

        # Compute the weighted energy values by including the population of each state
        weighted_energies = [
            energy * population for energy, population in zip(energies, populations)
        ]

        # Compute the entropy of the system
        boltzmann_entropy = -kB * np.sum(populations * np.log(populations))

        # Add a property entry to the property archive
        entry_name = method_el
        self.properties.add(entry_name)

        # Compute the helmotz free energy for the ensamble
        self.properties[method_el].helmotz_free_energy = (
            np.sum([weighted_energies]) - temperature * boltzmann_entropy,
        )


# class MDTrajectory:
#     """
#     Iterator class for MD trajectories. Data is computed only when accessing the
#     elements of the object (via __getitem__ or __iter__)

#     Parameters
#     ----------
#     traj_path : str
#         path (prefix) of the trajectory files used for creating the MD run
#     method : str
#         level of theory at which the simulation was ran
#     """
#     def __init__(self, traj_filepath: str, method: str) -> None:

#         self.name = traj_filepath
#         self.method = method

#         self.md_out = f"MD_data/{traj_filepath}_md.out"
#         self.geo_end = f"MD_data/{traj_filepath}_geo_end.xyz"

#         self.box_side = None
#         if os.path.exists(f"MD_data/{traj_filepath}.pbc"):
#             with open(f"MD_data/{traj_filepath}.pbc") as f:
#                 self.box_side = float(f.read())

#         with open(self.geo_end, "r") as f:

#             self.atomcount = int(f.readline())
#             first_iter = int(f.readline().split()[-1])

#             for line in f:

#                 if "iter" in line:
#                     second_iter = int(line.split()[-1])
#                     self.mdrestartfreq = second_iter - first_iter
#                     break

#             for line in reversed(list(f)):
#                 if "iter" in line:
#                     nframes = int(line.split()[-1])
#                     break

#         self.frames = list(range(0, (nframes // self.mdrestartfreq) + 1))

#     def __iter__(self):
#         for index in self.frames:
#             with open(self.geo_end, "r") as geo_end:
#                 yield self.__getitem__(index, geo_end)

#     def __getitem__(self, index, geo_end: TextIOWrapper = None):
#         """returns the System object corresponding to the requested index/MD step

#         Parameters
#         ----------
#         index : int
#             MD step in the simulation (as list index, not ACTUAL MD iter number)
#         geo_end : TextIOWrapper, optional
#             if __getitem__ is called from __iter__, allows passing the geo_end.xyz file so
#             it is open only once

#         Returns
#         -------
#         system
#             System object corresponding to the requested frame in the MD simulation
#         """

#         MDindex = self.frames[index] * self.mdrestartfreq
#         close_flag = False

#         if not geo_end:
#             close_flag = True
#             geo_end = open(self.geo_end, "r")

#         with open(f"{self.name}_{MDindex}.xyz", "w+") as out:

#             self.atomcount = int(geo_end.readline())
#             out.write(f"{self.atomcount}\n")
#             start = None

#             for i, line in enumerate(geo_end):

#                 if line == f"MD iter: {MDindex}\n":
#                     start = i
#                     out.write(line)

#                 if start and i > start and i < start + self.atomcount + 1:
#                     out.write(
#                         f"{line.split()[0]}\t{line.split()[1]}\t{line.split()[2]}\t{line.split()[3]}\n"
#                     )

#                 if start and i > start + self.atomcount + 1:
#                     break

#         if close_flag:
#             geo_end.close()

#         # !!! implement bytestream input for System !!!
#         system = System(f"{self.name}_{MDindex}.xyz", box_side=self.box_side)
#         os.remove(f"{self.name}_{MDindex}.xyz")

#         with open(self.md_out, "r") as md_out:
#             found = False
#             for line in md_out:
#                 if line == f"MD step: {MDindex}\n":
#                     found = True
#                 if "Total MD Energy" in line and found:
#                     system.energies[self.method] = Energies(
#                         method=self.method,
#                         electronic=line.split()[3],
#                         vibronic=None,
#                     )
#                     break

#         return system

#     def __len__(self):
#         return len(self.frames)
