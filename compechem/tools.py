import os
import shutil
import pickle
from rdkit import Chem
from compechem.systems import System
import logging

logger = logging.getLogger(__name__)


def generate_inchikey(molfile: str):
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


def generate_inchi(molfile: str):
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


def add_flag(mol: System, flag: str):
    """Adds a warning flag to a System object.

    Parameters
    ----------
    mol : System object
        System object to which the flag will be added
    flag : str
        String representing the warning which needs to be added
    """
    mol.flags.append(flag)
    return


def dump(obj, filename: bool = None):
    """Generates a pickle file containing an object 

    Parameters
    ----------
    obj : anything
        Object to dump to pickle file. Can be anything, including individual System objects
    filename : str
        string containing the filename of the pickle file.

    Returns
    -------
    Saves a pickle file containing the input object
    """

    if filename is None and type(obj) == System:
        filename = f"{obj.name}.pickle"

    pickle.dump(obj, open(filename, "wb"))


def save_ext(ext: list, output_dir: str):
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


def process_output(
    mol: System, method: str, calc: str, charge: int = None, spin: int = None,
):
    """Processes the output of a calculation, copying the output files to a safe directory in the
    parent directory tree, and cleans the temporary directory if requested.

    Parameters
    ----------
    mol : System object
        Systems processed in the calculation
    method : str
        level of theory for the calculation
    calc : str
        Type of calculation
    charge : int, optional
        Charge of the molecule in the calculation
    spin : int, optional
        Spin of the molecule in the calculation
    """

    if charge is None:
        charge = mol.charge
    if spin is None:
        spin = mol.spin

    os.makedirs("../output_files", exist_ok=True)
    if os.path.exists("output.out"):
        shutil.copy(
            "output.out", f"../output_files/{mol.name}_{charge}_{spin}_{method}_{calc}.out",
        )

    if os.path.exists("output.err"):
        os.makedirs("../error_files", exist_ok=True)
        shutil.copy(
            "output.err", f"../error_files/{mol.name}_{charge}_{spin}_{method}_{calc}.err",
        )


def cyclization_check(start_file: str, end_file: str):
    """Checks if a cyclization has occurred (e.g., during a
    geometry optimization), or if a ring opening has occurred.

    Parameters
    ----------
    start_file : str
        .xyz file containing the starting geometry
    end_file : str
        .xyz file containing the final geometry

    Returns
    -------
    If a change in the number of rings is detected, returns True.
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

    if start_rings != end_rings:
        return True
    else:
        return False


def dissociation_check():
    """Checks if a dissociation has occurred (e.g., during a
    geometry optimization)

    Requirements
    ----------
    A .mol file in the current directory.

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
        return True
    else:
        return False


def split_multixyz(
    mol: System, file: str, suffix: str, charge: int = None, spin: int = None
):
    """Splits a .xyz file containing multiple structures into individual structures.

    Parameters
    ----------
    mol : System object
        Input molecule, giving the charge/spin (if not defined) and name of the output molecules
    file : str
        .xyz file containing the multiple structures
    suffix : str
        suffix to add to the new molecule names
    charge : int, optional
        Charge of the output molecules, by default the same as the input molecule
    spin : int, optional
        Spin of the output molecules, by default the same as the input molecule

    Returns
    -------
    molecules_list : list
        List containing the individual System object, whose structure is taken from the .xyz file
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
            with open(f"{mol.name}_{suffix}{num}.xyz", "w") as out:
                for _ in range(molsize + 2):
                    out.write(line)
                    line = f.readline()
            molecules_list.append(System(f"{mol.name}_{suffix}{num}.xyz", charge, spin))
            num += 1

    return molecules_list


def save_dftb_trajectory(output_prefix):
    """Saves the geo_end.xyz and md.out files to a temporary directory where an MDTrajectory
    object can go read the data it needs.

    Parameters
    ----------
    output_prefix : str
        name of the output trajectory files prefix
    """

    os.makedirs("../MD_data", exist_ok=True)

    shutil.move("md.out", f"../MD_data/{output_prefix}_md.out")
    shutil.move("geo_end.xyz", f"../MD_data/{output_prefix}_geo_end.xyz")


def parse_dftb_trajectory(output_name):
    """Parses a geo_end.xyz trajectory and an md.out file to export a single trajectory
    file also containing the energies for all frames

    Parameters
    ----------
    output_name : str
        name of the output trajectory file
    """

    os.makedirs("../MD_trajectories", exist_ok=True)

    energies = []
    with open("md.out", "r") as f:
        for line in f:
            if "Total MD Energy" in line:
                energies.append(float(line.split()[3]))

    with open("geo_end.xyz", "r") as inp:
        with open(f"../MD_trajectories/{output_name}", "w") as out:
            for linenum, line in enumerate(inp):
                if linenum == 0:
                    atomcount = int(line)
                if linenum % (atomcount + 2) == 0:
                    out.write(line)
                if linenum % (atomcount + 2) == 1:
                    md_iter = line.rstrip("\n")
                    out.write(f"  {md_iter}\tEnergy: {energies.pop(0)} Eh\n")
                if linenum % (atomcount + 2) > 1:
                    out.write(
                        f"{line.split()[0]}\t{line.split()[1]}\t{line.split()[2]}\t{line.split()[3]}\n"
                    )

