import os
from rdkit import Chem
from tempfile import TemporaryDirectory

from rdkit.Chem import inchi

parent_dir = os.getcwd()
parent_dir = "/mnt/f/Downloads"


class Molecule:
    def __init__(self, charge=0, spin=1) -> None:
        self.moleculename: str
        self.charge: int = charge
        self.spin: int = spin
        self.atomnumber: int
        self.geometry: list = []
        # Energies
        self.total_energy: float = 0
        self.vibronic_contribution: float = 0
        self.spe_energy: float = 0


def generate_inchikey(molfile):
    mol = Chem.MolFromMolFile(molfile)
    inchikey = Chem.MolToInchiKey(mol)
    return inchikey


def read_xyz(xyz_file):

    moleculename = os.path.basename(xyz_file).strip(".xyz")

    geometry = []

    with open(xyz_file, "r") as file:
        for linenum, line in enumerate(file):
            if linenum == 0:
                atomnumber = line
            if linenum > 1 and linenum < int(atomnumber) + 2:
                geometry.append(line)

    return moleculename, atomnumber, geometry


def write_xyz(molecule, xyz_file):
    with open(xyz_file, "w") as file:
        file.write(str(molecule.atomnumber))
        file.write("\n")
        for line in molecule.geometry:
            file.write(line)


def update_geometry(xyz_file):

    geometry = []

    with open(xyz_file, "r") as file:
        for linenum, line in enumerate(file):
            if linenum == 0:
                atomnumber = line
            if linenum > 1 and linenum < int(atomnumber) + 2:
                geometry.append(line)

    return geometry


def cyclization_check(molecule, start_file, end_file):

    os.system(f"crest --testtopo {start_file} > start_topo.out 2>> start_topo.out")
    os.system(f"crest --testtopo {end_file} > end_topo.out 2>> end_topo.out")

    with open("start_topo.out", "r") as out:
        start_rings = 0
        for line in out:
            if "Total number of rings in the system" in line:
                start_rings = int(line.split()[-1])
                break

    with open("end_topo.out", "r") as out:
        end_rings = 0
        for line in out:
            if "Total number of rings in the system" in line:
                end_rings = int(line.split()[-1])
                break
            if "No (valid) input file!" in line:
                break

    if start_rings < end_rings:
        print(f"ERROR: Cyclization spotted for {molecule.moleculename}.\n")


def tautomer_search(molecule, nproc=1):

    print(f"INFO: {molecule.moleculename} - CREST tautomer search")

    with TemporaryDirectory(
        prefix=molecule.moleculename + "_", suffix="_crestTAUT", dir=os.getcwd()
    ) as wdir:

        os.chdir(wdir)
        write_xyz(molecule, "geom.xyz")

        os.system(
            f"crest geom.xyz --alpb water --charge {molecule.charge} --uhf {molecule.spin-1} --mquick --fstrict --tautomerize -T {nproc} > output.out 2>> output.out"
        )

        cyclization_check(molecule, "geom.xyz", "tautomers.xyz")

        molecule.geometry = update_geometry("tautomers.xyz")

        os.chdir(parent_dir)


def opt_xtb(molecule, conformers=True, nproc=1):

    print(f"INFO: {molecule.moleculename} - xTB OPT")

    with TemporaryDirectory(
        prefix=molecule.moleculename + "_", suffix="_xtbOPT", dir=os.getcwd()
    ) as wdir:

        os.chdir(wdir)
        write_xyz(molecule, "geom.xyz")

        os.system(
            f"xtb geom.xyz --alpb water --charge {molecule.charge} --uhf {molecule.spin-1} --ohess -P {nproc} > output.out 2>> output.out"
        )

        if conformers is True:
            os.system(
                f"crest xtbopt.xyz --alpb water --chrg {molecule.charge} --uhf {molecule.spin-1} --mquick -T {nproc} > conformers.out 2>> conformers.out"
            )
            os.system(
                f"xtb crest_best.xyz --alpb water --chrg {molecule.charge} --uhf {molecule.spin-1} --ohess -P {nproc} > output.out 2>> output.out"
            )

        # Checking for dissociation
        mol_file = [f for f in os.listdir(".") if f.endswith(".mol")][-1]
        end_mol = Chem.MolFromMolFile(
            mol_file, sanitize=False, removeHs=False, strictParsing=False
        )
        end_smiles = Chem.MolToSmiles(end_mol)
        if "." in end_smiles:
            print(
                f"ERROR: {molecule.moleculename} (charge {molecule.charge} spin {molecule.spin}) has undergone dissociation.\n"
            )

        cyclization_check(molecule, "geom.xyz", "xtbopt.xyz")

        with open("output.out", "r") as out:
            for line in out:
                if "total free energy" in line:
                    molecule.total_energy = float(line.split()[-3])
                if "G(RRHO) contrib." in line:
                    molecule.vibronic_contribution = float(line.split()[-3])

        molecule.geometry = update_geometry("xtbopt.xyz")

        os.chdir(parent_dir)


def spe_ccsd(molecule, nproc=1, maxcore=7500):

    print(f"INFO: {molecule.moleculename} - CCSD SPE")

    with TemporaryDirectory(
        prefix=molecule.moleculename + "_", suffix="_ccsdSPE", dir=os.getcwd()
    ) as wdir:

        os.chdir(wdir)
        write_xyz(molecule, "geom.xyz")

        with open("input.inp", "w") as inp:
            inp.write(
                f"%pal nproc {nproc} end\n"
                f"%maxcore {maxcore}\n"
                "! DLPNO-CCSD ano-pVTZ\n"
                "! RIJCOSX AutoAux\n"
                "%CPCM\n"
                "  SMD True\n"
                '  SMDsolvent "water"\n'
                "end\n"
                f"* xyzfile {molecule.charge} {molecule.spin} geom.xyz\n"
            )

        os.system("$ORCADIR/orca input.inp > output.out")

        with open("output.out", "r") as out:
            for line in out:
                if "FINAL SINGLE POINT ENERGY" in line:
                    molecule.spe_energy = float(line.split()[-3])
                    molecule.total_energy = (
                        molecule.spe_energy + molecule.vibronic_contribution
                    )

        os.chdir(parent_dir)


molecule = Molecule(0, 1)

molecule.moleculename, molecule.atomnumber, molecule.geometry = read_xyz(
    parent_dir + "/01_2,6-dimethoxyphenol.xyz"
)

for num, line in enumerate(molecule.geometry):
    print(num, line)

opt_xtb(molecule, nproc=8)

for num, line in enumerate(molecule.geometry):
    print(num, line)

tautomer_search(molecule, nproc=8)

for num, line in enumerate(molecule.geometry):
    print(num, line)
