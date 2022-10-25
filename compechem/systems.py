from io import TextIOWrapper
import os
import numpy as np
import logging
from itertools import chain

logger = logging.getLogger(__name__)


class Energies:
    """Molecular energies in Hartree"""

    def __init__(
        self,
        method: str = None,
        electronic: float = None,
        vibronic: float = None,
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

    def __str__(self):
        return f"method: {self.method}, el: {self.electronic} Eh, vib: {self.vibronic} Eh"


class Properties:
    """Class containing system properties (such as pKa)."""

    def __init__(self):
        self.pka: dict = {}


class System:
    """System object.

    Attributes
    ----------
    name : str
        name of the system, taken from the .xyz file
    charge : int
        total charge of the system
    spin : int
        total spin of the system (2S+1)
    atomcount : int
        number of atoms contained in the system
    geometry : list
        list containing the atomic coordinates of the system
    box_side : float, optional
        for periodic systems, defines the length (in Å) of the box side
    energies : dict
        dictionary containing the electronic/vibronic energies of the system,
        calculated at various levels of theory
    properties : dict
        dictionary containing the properties of the system, such as pKa
    flags : list
        list containing all "warning" flags which might be encountered during calculations.
    """

    def __init__(
        self,
        xyz_file: str,
        charge: int = 0,
        spin: int = 1,
        box_side: float = None,
    ) -> None:
        """
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
        """

        self.name = os.path.basename(xyz_file).strip(".xyz")
        self.charge: int = charge
        self.spin: int = spin

        self.atomcount: int = None
        self.geometry: list = []

        self.box_side = box_side

        if self.box_side:
            self.periodic = True
        else:
            self.periodic = False

        self.velocities: list = []

        self.flags: list = []

        self.energies: dict = {}
        self.properties: Properties = Properties()

        self.update_geometry(xyz_file)

    def __str__(self):
        info = f"=== System: {self.name} === \n"
        info += f"\nNumber of atoms: {self.atomcount}\n"
        info += f"Charge: {self.charge}\n"
        info += f"Spin: {self.spin}\n"
        info += "\n--- Warnings ---\n"
        for warning in self.flags:
            info += f"{warning}\n"
        info += "\n--- Energies (Eh) --- \n"
        for method in self.energies:
            info += f"\n* Method: {method}\n"
            info += f"Electronic: {self.energies[method].electronic} Eh\n"
            info += f"Vibronic: {self.energies[method].vibronic} Eh\n"
        info += "\n--- Coordinates (Å) --- \n\n"
        for line in self.geometry:
            info += f"{line}"
        info += "\n--- Velocities (Å/ps) --- \n\n"
        for line in self.velocities:
            info += f"{line}"

        return info

    def write_xyz(self, xyz_file: str):
        """Writes the current geometry to a .xyz file.

        Parameters
        ----------
        xyz_file : str
            path to the output .xyz file
        """
        with open(xyz_file, "w") as file:
            file.write(str(self.atomcount))
            file.write("\n\n")
            for line in self.geometry:
                file.write(line)

    def write_gen(self, gen_file: str, box_side: float = None):
        """Writes the current geometry to a .gen file.

        Parameters
        ----------
        gen_file : str
            path to the output .gen file
        box_side : float, optional
            for periodic systems, defines the length (in Å) of the box side
        """

        if box_side is None:
            box_side = self.box_side

        with open(gen_file, "w") as file:
            if self.periodic:
                file.write(f" {str(self.atomcount)} S\n")
            else:
                file.write(f" {str(self.atomcount)} C\n")
            atom_types = []
            for line in self.geometry:
                if line.split()[0] not in atom_types:
                    atom_types.append(line.split()[0])
            for atom in atom_types:
                file.write(f" {atom}")
            file.write("\n")
            i = 1
            for line in self.geometry:
                for index, atom in enumerate(atom_types):
                    if line.split()[0] == atom:
                        file.write(f"{i} {line.replace(atom, str(index + 1))}")
                        i += 1
            if self.periodic:
                file.write(f" 0.000 0.000 0.000\n")
                file.write(f" {box_side} 0.000 0.000\n")
                file.write(f" 0.000 {box_side} 0.000\n")
                file.write(f" 0.000 0.000 {box_side}")

    def update_geometry(self, xyz_file: str):
        """Updates the current geometry from an external .xyz file

        Parameters
        ----------
        xyz_file : str
            path with the .xyz file of the geometry containing the
            new coordinates
        """
        self.geometry = []

        with open(xyz_file, "r") as f:
            numlines = sum(1 for line in f)

        with open(xyz_file, "r") as f:
            for linenum, line in enumerate(f):
                if linenum == 0:
                    self.atomcount = int(line)
                last_geom_line = numlines - (self.atomcount + 2)
                if linenum > last_geom_line + 1 and linenum < numlines and len(line) != 0:
                    self.geometry.append(
                        f"{line.split()[0]}\t{line.split()[1]}\t{line.split()[2]}\t{line.split()[3]}\n"
                    )
                    if len(line.split()) > 4:
                        self.velocities.append(
                            f"{line.split()[0]}\t{line.split()[-3]}\t{line.split()[-2]}\t{line.split()[-1]}\n"
                        )


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
    container : list
        iterable returning System objects generator for creating the Ensemble.
    """

    def __init__(self, systems_list) -> None:
        """
        Parameters
        ----------
        systems_list : any iterable returning System objects
            generator for creating the Ensemble.
            Options:
                - list of System objects
                - MDTrajectory iterator
        """

        self.name = systems_list[0].name
        self.atomcount = systems_list[0].atomcount

        self.container = [systems_list]

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
        """append more Systems to the ensemble from an iterator object (e.g., MDTrajectory)

        Parameters
        ----------
        iterator : iterable
            iterator object (e.g., MDTrajectory object)
        """
        self.container.append(iterator)

    # NOT WORKING - FIX IT!
    def read_energies(self, method):
        """reads energies from trajectory file (parsed by tools.save_dftb_trajectory()) and
        calculates the average energy for the given trajectory

        CURRENTLY NOT WORKING!!!

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
    def __init__(self, traj_filepath, method):
        """Iterator class for MD trajectories. Data is computed only when accessing the
        elements of the object (via __getitem__ or __iter__)

        Parameters
        ----------
        traj_path : str
            path (prefix) of the trajectory files used for creating the MD run
        method : str
            level of theory at which the simulation was ran
        """
        self.name = traj_filepath
        self.method = method

        self.md_out = f"MD_data/{traj_filepath}_md.out"
        self.geo_end = f"MD_data/{traj_filepath}_geo_end.xyz"

        self.periodic = False
        self.box_side = None
        if os.path.exists(f"MD_data/{traj_filepath}.pbc"):
            self.periodic = True
            with open(f"MD_data/{traj_filepath}.pbc") as f:
                self.box_side = float(f.read())

        with open(self.geo_end, "r") as f:

            self.atomcount = int(f.readline())
            first_iter = int(f.readline().split()[-1])

            for line in f:

                if "iter" in line:
                    second_iter = int(line.split()[-1])
                    self.mdrestartfreq = second_iter - first_iter
                    break

            for line in reversed(list(f)):
                if "iter" in line:
                    nframes = int(line.split()[-1])
                    break

        self.frames = list(range(0, (nframes // self.mdrestartfreq) + 1))

    def __iter__(self):
        for index in self.frames:
            with open(self.geo_end, "r") as geo_end:
                yield self.__getitem__(index, geo_end)

    def __getitem__(self, index, geo_end: TextIOWrapper = None):
        """returns the System object corresponding to the requested index/MD step

        Parameters
        ----------
        index : int
            MD step in the simulation (as list index, not ACTUAL MD iter number)
        geo_end : TextIOWrapper, optional
            if __getitem__ is called from __iter__, allows passing the geo_end.xyz file so
            it is open only once

        Returns
        -------
        system
            System object corresponding to the requested frame in the MD simulation
        """

        MDindex = self.frames[index] * self.mdrestartfreq
        close_flag = False

        if not geo_end:
            close_flag = True
            geo_end = open(self.geo_end, "r")

        with open(f"{self.name}_{MDindex}.xyz", "w+") as out:

            self.atomcount = int(geo_end.readline())
            out.write(f"{self.atomcount}\n")
            start = None

            for i, line in enumerate(geo_end):

                if line == f"MD iter: {MDindex}\n":
                    start = i
                    out.write(line)

                if start and i > start and i < start + self.atomcount + 1:
                    out.write(
                        f"{line.split()[0]}\t{line.split()[1]}\t{line.split()[2]}\t{line.split()[3]}\n"
                    )

                if start and i > start + self.atomcount + 1:
                    break

        if close_flag:
            geo_end.close()

        # !!! implement bytestream input for System !!!
        system = System(
            f"{self.name}_{MDindex}.xyz", periodic=self.periodic, box_side=self.box_side
        )
        os.remove(f"{self.name}_{MDindex}.xyz")

        with open(self.md_out, "r") as md_out:
            found = False
            for line in md_out:
                if line == f"MD step: {MDindex}\n":
                    found = True
                if "Total MD Energy" in line and found:
                    system.energies[self.method] = Energies(
                        method=self.method,
                        electronic=line.split()[3],
                        vibronic=None,
                    )
                    break

        return system

    def __len__(self):
        return len(self.frames)
