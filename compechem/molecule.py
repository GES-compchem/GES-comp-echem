import os


class Molecule:
    def __init__(self, xyz_file, charge=0, spin=1) -> None:
        self.name: str
        self.charge: int = charge
        self.spin: int = spin
        self.atomcount: int
        self.geometry: list = []

        # Energies
        self.energies = {
            "total_energy": 0.0,
            "vibronic_contribution": 0.0,
            "spe_energy": 0.0,
        }

        self.name = os.path.basename(xyz_file).strip(".xyz")

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    self.atomcount = int(line)
                if linenum > 1 and linenum < self.atomcount + 2:
                    self.geometry.append(line)

    def update_geometry(self, xyz_file):

        self.geometry = []

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    atomcount = int(line)
                if linenum > 1 and linenum < atomcount + 2:
                    self.geometry.append(line)

