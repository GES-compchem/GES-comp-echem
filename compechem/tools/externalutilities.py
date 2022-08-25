import os
from compechem.systems import System
import logging

logger = logging.getLogger(__name__)


def split_multixyz(
    mol: System, file: str, suffix: str = "", charge: int = None, spin: int = None
):
    """Splits a .xyz file containing multiple structures into individual structures.

    Parameters
    ----------
    mol : System object
        Input molecule, giving the charge/spin (if not defined) and name of the output molecules
    file : str
        .xyz file containing the multiple structures
    suffix : str, optional
        suffix to add to the new molecule names. By default, empty.
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


def compress_dftb_trajectory(filename, md_out="md.out", geo_xyz="geo_end.xyz"):
    """Parses a geo_end.xyz trajectory and an md.out file to export a single compressed
    trajectory file also containing the energies for all frames

    Parameters
    ----------
    filename : str
        name of the output trajectory files
    md_out : str, optional
        path to the md.out file containing energy info (by default, ./md.out)
    geo_xyz : str, optional
        path to the geo_end.xyz file containing energy info (by default, ./geo_end.xyz)
    """

    logger.info(f"Parsing trajectory file: {geo_xyz}")
    logger.info(f"Parsing energy file: {md_out}")

    energies = []
    with open(md_out, "r") as f:
        for line in f:
            if "Total MD Energy" in line:
                energies.append(float(line.split()[3]))

    with open(geo_xyz, "r") as inp:
        with open(f"{filename}.xyz", "w") as out:
            for linenum, line in enumerate(inp):
                if linenum == 0:
                    atomcount = int(line)
                if linenum % (atomcount + 2) == 0:
                    out.write(f"{line.split()[0]}\n")
                if linenum % (atomcount + 2) == 1:
                    out.write(f"Step: {line.split()[0]} Energy: {energies.pop(0)} Eh\n")
                if linenum % (atomcount + 2) > 1:
                    out.write(
                        f"{line.split()[0]} {round(float(line.split()[1]),3)} {round(float(line.split()[2]),3)} {round(float(line.split()[3]),3)}\n"
                    )
    logger.info(f"Compressing MD trajectory to {filename}.zip")
    os.system(f"zip {filename}.zip {filename}.xyz")
