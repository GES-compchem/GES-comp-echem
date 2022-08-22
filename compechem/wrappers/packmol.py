import os, shutil, sh
import logging
from tempfile import mkdtemp
from compechem.systems import System
from compechem.tools import process_output

logger = logging.getLogger(__name__)


def packmol_cube(
    solute: System,
    solvent: System,
    nsolv: int = None,
    target_dens: float = None,
    cube_side: float = None,
):
    """Interface for running packmol to generate a cubic solvation box.

    Parameters
    ----------
    solute : str
        System object of the solute molecule
    solvent : str
        System object of the solvent molecule

    # two of the following three parameters are required, the third will be calculated
    # based on the other two
    nsolv : int, optional
        number of solvent molecules to add in the box
    target_dens : float, optional
        target density for the solvated box, in g/L
    cube_side : float, optional
        length of the box side, in Å

    Returns
    -------
    Returns a System object
    """

    avogadro = 6.0221408e23

    mol_weights = {
        "H": 1.00797,
        "B": 10.81,
        "C": 12.011,
        "N": 14.0067,
        "O": 15.9994,
        "F": 18.998403,
        "S": 32.06,
        "Cl": 35.453,
        "Br": 79.904,
        "I": 126.9045,
    }

    if nsolv and target_dens and cube_side:
        logger.error(
            "At least one of ( nsolv | target_dens | cube_side ) must be left out."
        )
        return None

    elif nsolv and target_dens:

        solvent_grams = (
            sum([mol_weights[atom[0]] for atom in solvent.geometry]) * nsolv / avogadro
        )  # g
        solute_grams = (
            sum([mol_weights[atom[0]] for atom in solute.geometry]) / avogadro
        )  # g

        target_volume = (solvent_grams + solute_grams) / target_dens  # L

        cube_side = (target_volume / 1e-27) ** (1.0 / 3)

    elif nsolv and cube_side:

        solvent_grams = (
            sum([mol_weights[atom[0]] for atom in solvent.geometry]) * nsolv / avogadro
        )  # g
        solute_grams = (
            sum([mol_weights[atom[0]] for atom in solute.geometry]) / avogadro
        )  # g

        volume = (cube_side**3) * 1e-27  # L

        target_dens = (solvent_grams + solute_grams) / volume

    elif target_dens and cube_side:

        volume = (cube_side**3) * 1e-27  # L

        solute_grams = (
            sum([mol_weights[atom[0]] for atom in solute.geometry]) / avogadro
        )  # g

        target_solv_weight = (target_dens * volume) - solute_grams  # g

        nsolv = round(
            target_solv_weight
            * avogadro
            / sum([mol_weights[atom[0]] for atom in solvent.geometry])
        )  # n° of molecules

        solvent_grams = (
            sum([mol_weights[atom[0]] for atom in solvent.geometry]) * nsolv / avogadro
        )  # g

        # recalculating density with actual number of solvent molecules
        target_dens = (solvent_grams + solute_grams) / volume

    else:
        logger.error(
            "At least two of ( nsolv | target_dens | cube_side ) must be provided."
        )
        return None

    logger.info(f"{solute.name} - Generating solvation box with {nsolv} {solvent.name}s")
    # Also print EFFECTIVE DENSITY!
    logger.debug(
        f"Packmol solvated {solute.name} - cubic box with {nsolv} {solvent.name} molecules, side {cube_side} Å, density {target_dens} g/L."
    )

    tdir = mkdtemp(
        prefix=f"{solute.name}_{nsolv}{solvent.name}s" + "_",
        suffix=f"_packmol",
        dir=os.getcwd(),
    )

    with sh.pushd(tdir):

        solute.write_xyz(f"{solute.name}.xyz")
        solvent.write_xyz(f"{solvent.name}.xyz")

        with open("input.inp", "w") as f:
            f.write(
                "tolerance 2.0\n"
                f"filetype xyz\n\n"
                f"output {solute.name}_{nsolv}{solvent.name}s.xyz\n\n"
                f"structure {solvent.name}.xyz\n"
                f"  number {nsolv}\n"
                "  resnumbers 3\n"
                f"  inside cube 0. 0. 0. {cube_side-2}\n"
                "end structure\n\n"
                f"structure {solute.name}.xyz\n"
                "  number 1\n"
                "  resnumbers 3\n"
                "  center\n"
                f"  fixed {(cube_side-2)/2} {(cube_side-2)/2} {(cube_side-2)/2} 0. 0. 0.\n"
                "end structure\n"
            )

        os.system("packmol < input.inp > output.out")

        solvated_molecule = System(
            f"{solute.name}_{nsolv}{solvent.name}s.xyz",
            periodic=True,
            box_side=cube_side,
        )

        process_output(mol=solute, method="packmol", calc=f"{nsolv}{solvent.name}_cube")
        shutil.rmtree(tdir)

    return solvated_molecule
