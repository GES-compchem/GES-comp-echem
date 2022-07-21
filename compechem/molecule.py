import os


class Energies:
    """Molecular energies in Hartree
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

    def __str__(self):
        return f"method: {self.method}, el: {self.electronic} Eh, vib: {self.vibronic} Eh"


class Properties:
    """Class containing molecule properties (such as pKa)."""

    def __init__(self):
        self.pka: dict = {}


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
    flags : list
        list containing all "warning" flags which might be encountered during calculations.
    """

    def __init__(
        self,
        xyz_file: str,
        charge: int = 0,
        spin: int = 1,
        geom_type: str = "C",
        box_side: float = None,
    ) -> None:
        """
        Parameters
        ----------
        xyz_file : str
            path with the .xyz file containing the molecule geometry
        charge : int, optional
            total charge of the molecule. Defaults to 0 (neutral)
        spin : int, optional
            total spin of the molecule. Defaults to 1 (singlet)
        geom_type : str, optional
            type of geometry for the molecule. 
                Options: 
                - C (default) = cluster (single molecule)
                - S = Supercell (periodic system)
        box_side : float, optional
            for periodic systems, defines the length (in Å) of the box side
        """

        self.name = os.path.basename(xyz_file).strip(".xyz")
        self.charge: int = charge
        self.spin: int = spin

        self.atomcount: int = None
        self.geometry: list = []

        self.geom_type = geom_type
        self.box_side = box_side

        self.flags: list = []

        self.energies: dict = {}
        self.properties: Properties = Properties()

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    self.atomcount = int(line)
                if linenum > 1 and linenum < self.atomcount + 2:
                    self.geometry.append(line)

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
            file.write(f" {str(self.atomcount)} {self.geom_type}\n")
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
            if self.geom_type == "S":
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

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    self.atomcount = int(line)
                if linenum > 1 and linenum < self.atomcount + 2:
                    self.geometry.append(line)
