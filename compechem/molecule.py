import os


class Molecule:
    """Molecule object.

    Attributes
    ----------
    name : str
        name of the molecule, taken from the .xyz file
    charge : int
        total charge of the molecule
    spin : int
        total spin of the molecule (2S+1)
    atomcount : int
        number of atoms contained in the molecule
    geometry : list
        list containing the atomic coordinates of the molecule
    energies : dict
        dictionary containing the electronic/vibronic energies of the molecule,
        calculated at various levels of theory
    """

    def __init__(self, xyz_file, charge=0, spin=1) -> None:
        """
        Parameters
        ----------
        xyz_file : str
            path with the .xyz file containing the molecule geometry
        charge : int, optional
            total charge of the molecule. Defaults to 0 (neutral)
        spin : int, optional
            total spin of the molecule. Defaults to 1 (singlet)
        """

        self.name = os.path.basename(xyz_file).strip(".xyz")
        self.charge: int = charge
        self.spin: int = spin
        self.atomcount: int = None
        self.geometry: list = []

        self.energies: dict = {}

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    self.atomcount = int(line)
                if linenum > 1 and linenum < self.atomcount + 2:
                    self.geometry.append(line)

    class Energies:
        """Molecular energies.
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

    def write_xyz(self, xyz_file):
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

    def update_geometry(self, xyz_file):
        """Updates the current geometry from an external .xyz file

        Parameters
        ----------
        xyz_file : str
            path with the .xyz file of the geometry containing the 
            new coordinates
        """
        self.geometry = []

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    self.atomcount = int(line)
                if linenum > 1 and linenum < self.atomcount + 2:
                    self.geometry.append(line)
