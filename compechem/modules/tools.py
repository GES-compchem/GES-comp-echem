import os
import shutil
import pickle
from rdkit import Chem
from compechem.modules.orca import OrcaInput
from compechem.modules.xtb import XtbInput
from compechem.molecule import Molecule


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


def info(mol, print_geometry=True):
    """Prints information about the molecule

    Parameters
    ----------
    mol : Molecule object
        Molecule object
    print_geometry : bool
        prints atom coordinates, by default True

    Returns
    -------
    Prints to screen a summary with all the informations about the molecule.
    """
    print(f"\n === Molecule: {mol.name} === ")
    print(f"\nNumber of atoms: {mol.atomcount}")
    print(f"Charge: {mol.charge}")
    print(f"Spin: {mol.spin}")
    print("\n --- Energies (Eh) --- \n")
    for method in mol.energies:
        print(f"* Method: {method}")
        print(f"Electronic: {mol.energies[method].electronic}")
        print(f"Vibronic: {mol.energies[method].vibronic}")
    if print_geometry is True:
        print("\n --- Coordinates (Angstrom) --- ")
        for line in mol.geometry:
            print(line, end="")


def dump(obj, filename=None):
    """Generates a pickle file containing an object 

    Parameters
    ----------
    obj : anything
        Object to dump to pickle file. Can be anything, including individual Molecule objects
    filename : str
        string containing the filename of the pickle file.

    Returns
    -------
    Saves a pickle file containing the input object
    """

    if filename is None and type(obj) == Molecule:
        filename = f"{obj.name}.pickle"

    pickle.dump(obj, open(filename, "wb"))


def save_ext(ext, output_dir):
    """Saves all files matching a certain set of extensions

    Parameters
    ----------
    ext : list
        List containing the extensions of the files to be saved (in string format)
    output_dir : str
        Directory in which the files are to be saved

    Returns
    -------
    Saves all the files matching the given extensions to the output directory
    """

    os.makedirs(output_dir, exist_ok=True)

    for file in os.listdir("."):
        if os.path.splitext(file)[1] in ext:
            shutil.copy(file, output_dir + "/" + file)


def process_output(mol, calc, runtype, tdir, remove_tdir, parent_dir):
    """Processes the output of a calculation, copying the output files to a safe directory in the
    parent directory tree, and cleans the temporary directory if requested.

    Parameters
    ----------
    mol : Molecule object
        Molecules processed in the calculation
    calc : Input object or str
        Input (xtb/orca) object containing information for the calculation, or generic string
        if not one of the above.
    runtype : str
        Type of calculation
    tdir : str
        Temporary directory
    remove_tdir : bool
        If true, removes the temporary directory
    parent_dir : str
        Parent directory to return to after the calculation is done    

    """

    os.makedirs("../output_files", exist_ok=True)
    if calc is OrcaInput or XtbInput:
        shutil.copy(
            "output.out",
            f"../output_files/{mol.name}_{calc.charge}_{calc.spin}_{calc.method}_{runtype}.out",
        )
    else:
        shutil.copy(
            "output.out", f"../output_files/{mol.name}_{calc}_{runtype}.out",
        )

    if os.path.exists("output.err"):
        os.makedirs("../error_files", exist_ok=True)
        if calc is OrcaInput or XtbInput:
            shutil.copy(
                "output.err",
                f"../error_files/{mol.name}_{calc.charge}_{calc.spin}_{calc.method}_{runtype}.err",
            )
        else:
            shutil.copy(
                "output.err", f"../error_files/{mol.name}_{calc}_{runtype}.err",
            )

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)


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
    mol : Molecule object
        Molecule for which to check if a dissociation has happened.
        Requires a .mol file in the current directory.

    Returns
    -------
    True, if a dissociation is detected
    False, if a dissociation is not detected
    """

    mol_file = [f for f in os.listdir(".") if f.endswith(".mol")][-1]
    end_mol = Chem.MolFromMolFile(mol_file, sanitize=False, removeHs=False, strictParsing=False)
    end_smiles = Chem.MolToSmiles(end_mol)
    if "." in end_smiles:
        print(
            f"ERROR: {mol.name} (charge {mol.charge} spin {mol.spin}) has undergone dissociation.\n"
        )
        return True

    else:
        return False


def split_multixyz(mol, file, charge=None, spin=None):
    """Splits a .xyz file containing multiple structures into individual structures.

    Parameters
    ----------
    mol : Molecule object
        Input molecule, giving the charge/spin (if not defined) and name of the output molecules
    file : str
        .xyz file containing the multiple structures
    charge : int, optional
        Charge of the output molecules, by default the same as the input molecule
    spin : int, optional
        Spin of the output molecules, by default the same as the input molecule

    Returns
    -------
    molecules_list : list
        List containing the individual Molecule object, whose structure is taken from the .xyz file
    """

    if charge is None:
        charge = mol.charge
    if spin is None:
        spin = mol.spin

    with open(file, "r") as f:
        molsize = int(f.readline())

    molecules_list = []

    num = 1
    with open(file, "r") as f:
        line = f.readline()
        while line:
            with open(f"{mol.name}_{os.path.basename(file).strip('.xyz')}_{num}.xyz", "w") as out:
                for _ in range(molsize + 2):
                    out.write(line)
                    line = f.readline()
            molecules_list.append(
                Molecule(
                    f"{mol.name}_{os.path.basename(file).strip('.xyz')}_{num}.xyz", charge, spin
                )
            )
            num += 1

    return molecules_list
