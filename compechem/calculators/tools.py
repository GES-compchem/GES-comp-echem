import os
from rdkit import Chem


def generate_inchikey(molfile):
    mol = Chem.MolFromMolFile(molfile)
    inchikey = Chem.MolToInchiKey(mol)
    return inchikey


def info(mol):
    print(f"\nMolecule: {mol.name}")
    print(f"\tNumber of atoms: {mol.atomcount}")
    print(f"\tCharge: {mol.charge}")
    print(f"\tSpin: {mol.spin}")
    print("\nEnergies:")
    print(f"\t Total: {mol.energies['total']} Eh")
    print(f"\t Vibronic: {mol.energies['vibronic']} Eh")
    print(f"\t Electronic: {mol.energies['electronic']} Eh")
    print("\nCoordinates (Angstrom):")
    for line in mol.geometry:
        print(line, end="")
    print("\n")


def write_xyz(mol, xyz_file):
    with open(xyz_file, "w") as file:
        file.write(str(mol.atomcount))
        file.write("\n\n")
        for line in mol.geometry:
            file.write(line)


def cyclization_check(mol, start_file, end_file):

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
        print(f"ERROR: Cyclization spotted for {mol.name}.\n")


def dissociation_check(mol):

    mol_file = [f for f in os.listdir(".") if f.endswith(".mol")][-1]
    end_mol = Chem.MolFromMolFile(
        mol_file, sanitize=False, removeHs=False, strictParsing=False
    )
    end_smiles = Chem.MolToSmiles(end_mol)
    if "." in end_smiles:
        print(
            f"ERROR: {mol.name} (charge {mol.charge} spin {mol.spin}) has undergone dissociation.\n"
        )
        return True

    else:
        return False
