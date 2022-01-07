import os, shutil
from rdkit import Chem
from tempfile import mkdtemp

parent_dir = os.getcwd()


class Molecule:
    def __init__(self, xyz_file, charge=0, spin=1) -> None:
        self.moleculename: str
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
        self.moleculename = os.path.basename(xyz_file).strip(".xyz")

        self.geometry: list = []
        self.deprotonation_state: int = 0  # how many times it was deprotonated

        with open(xyz_file, "r") as file:
            for linenum, line in enumerate(file):
                if linenum == 0:
                    self.atomcount = line
                if linenum > 1 and linenum < int(self.atomcount) + 2:
                    self.geometry.append(line)


def generate_inchikey(molfile):
    mol = Chem.MolFromMolFile(molfile)
    inchikey = Chem.MolToInchiKey(mol)
    return inchikey


def write_xyz(molecule, xyz_file):
    with open(xyz_file, "w") as file:
        file.write(str(molecule.atomcount))
        file.write("\n")
        for line in molecule.geometry:
            file.write(line)


def update_geometry(molecule, xyz_file):

    molecule.geometry = []

    with open(xyz_file, "r") as file:
        for linenum, line in enumerate(file):
            if linenum == 0:
                atomcount = line
            if linenum > 1 and linenum < int(atomcount) + 2:
                molecule.geometry.append(line)


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

    wdir = mkdtemp(
        prefix=molecule.moleculename + "_", suffix="_crestTAUT", dir=os.getcwd()
    )

    os.chdir(wdir)
    write_xyz(molecule, "geom.xyz")

    os.system(
        f"crest geom.xyz --alpb water --charge {molecule.charge} --uhf {molecule.spin-1} --mquick --fstrict --tautomerize -T {nproc} > output.out 2>> output.out"
    )

    cyclization_check(molecule, "geom.xyz", "tautomers.xyz")

    update_geometry(molecule, "tautomers.xyz")

    shutil.rmtree(wdir)
    os.chdir(parent_dir)


def opt_xtb(molecule, conformers=True, nproc=1):

    print(f"INFO: {molecule.moleculename} - xTB OPT")

    wdir = mkdtemp(
        prefix=molecule.moleculename + "_", suffix="_xtbOPT", dir=os.getcwd()
    )

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
                molecule.energies["total_energy"] = float(line.split()[-3])
            if "G(RRHO) contrib." in line:
                molecule.energies["vibronic_contribution"] = float(line.split()[-3])

    update_geometry(molecule, "xtbopt.xyz")

    shutil.rmtree(wdir)
    os.chdir(parent_dir)


def spe_ccsd(molecule, nproc=1, maxcore=7500):

    print(f"INFO: {molecule.moleculename} - CCSD SPE")

    wdir = mkdtemp(
        prefix=molecule.moleculename + "_", suffix="_ccsdSPE", dir=os.getcwd()
    )

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
                molecule.energies["spe_energy"] = float(line.split()[-3])
                molecule.energies["total_energy"] = (
                    molecule.energies["spe_energy"]
                    + molecule.energies["vibronic_contribution"]
                )

    shutil.rmtree(wdir)
    os.chdir(parent_dir)


def spe_b97(molecule, nproc=1, maxcore=350):

    print(f"INFO: {molecule.moleculename} - B97 SPE")

    wdir = mkdtemp(
        prefix=molecule.moleculename + "_", suffix="_b97SPE", dir=os.getcwd()
    )

    os.chdir(wdir)
    write_xyz(molecule, "geom.xyz")

    with open("input.inp", "w") as inp:
        inp.write(
            f"%pal nproc {nproc} end\n"
            f"%maxcore {maxcore}\n"
            "! B97-D3 D3BJ def2-TZVP\n"
            "! RIJCOSX def2/J\n"
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
                molecule.energies["spe_energy"] = float(line.split()[-3])
                molecule.energies["total_energy"] = (
                    molecule.energies["spe_energy"]
                    + molecule.energies["vibronic_contribution"]
                )

    shutil.rmtree(wdir)
    os.chdir(parent_dir)


def deprotonate(molecule, nproc=1):

    print(f"INFO: {molecule.moleculename} - CREST deprotonation")

    wdir = mkdtemp(
        prefix=molecule.moleculename + "_", suffix="_crestDEPROT", dir=os.getcwd()
    )

    os.chdir(wdir)
    write_xyz(molecule, "geom.xyz")

    os.system(
        f"crest geom.xyz --alpb water --charge {molecule.charge} --uhf {molecule.spin-1} --deprotonate -T {nproc} > output.out 2>> output.out"
    )

    molsize = int(molecule.atomcount) - 1
    with open("deprotonated.xyz", "r") as f:
        numlines = int(sum(1 for line in f))

    j = 0
    deprotomer = 1
    with open("deprotonated.xyz", "r") as f:
        while j < numlines:
            geometry = []
            i = 0
            while i < molsize + 2:
                geometry.append(f.readline())
                i += 1
                j += 1
            molecule.deprotomers.append(geometry)
            deprotomer += 1

    shutil.rmtree(wdir)
    os.chdir(parent_dir)


def calculate_pka(molecule, nproc=1):

    print(
        f"INFO: calculating pKa for {molecule.moleculename}, charge {molecule.charge} spin {molecule.spin}"
    )

