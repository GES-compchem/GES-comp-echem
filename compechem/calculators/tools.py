import os
from rdkit import Chem


def generate_inchikey(molfile):
    """Generates InchiKey from a Mol file.

    Parameters
    ----------
    molfile : str
        path containing the Mol file

    Returns
    -------
    InchiKey
        InchiKey generated from the Mol file
    """
    mol = Chem.MolFromMolFile(molfile)
    inchikey = Chem.MolToInchiKey(mol)
    return inchikey


def generate_inchi(molfile):
    """Generates Inchi from a Mol file.

    Parameters
    ----------
    molfile : str
        path containing the Mol file

    Returns
    -------
    Inchi
        Inchi generated from the Mol file
    """
    mol = Chem.MolFromMolFile(molfile)
    inchi = Chem.MolToInchi(mol)
    return inchi


def info(mol):
    """Prints information about the molecule

    Parameters
    ----------
    mol : Molecule object
        Molecule object

    Returns
    -------
    Prints to screen a summary with all the informations about the molecule.
    """
    print(f"\nMolecule: {mol.name}")
    print(f"\tNumber of atoms: {mol.atomcount}")
    print(f"\tCharge: {mol.charge}")
    print(f"\tSpin: {mol.spin}")
    print("\nEnergies (Eh):")
    for method in mol.energies:
        print(f"\n\t Method: {method}")
        print(f"\t Electronic: {mol.energies[method].electronic}")
        print(f"\t Vibronic: {mol.energies[method].vibronic}")
    print("\nCoordinates (Angstrom):")
    for line in mol.geometry:
        print(line, end="")
    print("\n")


def cyclization_check(mol, start_file, end_file):
    """Checks if a cyclization has occurred (e.g., during a
    geometry optimization)

    Parameters
    ----------
    mol : Molecule object
        Molecule object
    start_file : str
        .xyz file containing the starting geometry
    end_file : str
        .xyz file containing the final geometry

    Returns
    -------
    If a cyclization is detected, prints an error message.
    """

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
    """Checks if a dissociation has occurred (e.g., during a
    geometry optimization)

    Parameters
    ----------
    mol : str
        Mol file containing the molecular structure

    Returns
    -------
    True, if a dissociation is detected
    False, if a dissociation is not detected
    """

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
