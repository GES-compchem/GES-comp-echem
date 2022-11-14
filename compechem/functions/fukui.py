from os import listdir
from os.path import join
from copy import deepcopy
from typing import List, Dict

from compechem.systems import System
from compechem.wrappers.orca import OrcaInput
from compechem.tools.cubetools import Cube


def compute_fukui(
    xyz_file: str,
    charge: int = 0,
    spins_states: List[int] = [2, 1, 2],
    method: str = "B3LYP",
    basis_set: str = "def2-TZVP",
) -> None:
    """
    Computes the Fukui f+, f- and f0 functions starting from a given input molecule. The functions
    are saved in Gaussian cube format and stored in the `output_densities` folder.    

    Parameters
    ----------
    xyz_file: str
        the path to the .xyz file encoding the geometry of the starting molecule.
    charge: int
        the charge associated to the molecule
    spin_states: List[int]
        the list of spin multeplicities associated with the molecule with one
        added electron, with the molecule as it is and with the molecule with one
        of its electrons removed.
    method: str
        the level of theory at which the calculations should be performed
    basis_set: str
        the basis-set to be used in the computation
    """
    #Define the molecule for which the Fukui functions must be computed 
    origin = System(xyz_file, charge=charge, spin=spins_states[1])

    # Define an Orca interface using the given method inputs
    orca = OrcaInput(method=method, basis_set=basis_set, aux_basis="")

    # Optimize the molecule
    orca.opt(origin, elden_cube=True, inplace=True)

    # Compute a single point for the molecules with the addition of one electron.
    cation = deepcopy(origin)
    cation.charge += 1
    cation.spin = spins_states[2]
    orca.spe(cation, elden_cube=True, inplace=True)

    # Compute a single point for the molecules with the subtraction of one electron.
    anion = deepcopy(origin)
    anion.charge -= 1
    anion.spin = spins_states[0]
    orca.spe(anion, elden_cube=True, inplace=True)

    # Load cubes from the output_densities folder
    cubes: Dict[int, Cube] = {}
    for file in listdir("./output_densities"):
        if file.endswith("eldens.cube"):
            cube = Cube.from_file(join("./output_densities", file))
            current_charge = int(file.split("_")[1])
            cubes[current_charge] = cube

    # Check if all the densities have been loaded correctly
    if len(cubes) != 3:
        raise RuntimeError(f"Three cube files expected, {len(cubes)} found.")
    
    # Compute the f+ Fukui function
    f_plus = cubes[charge+1] - cubes[charge]
    f_plus.save(join("./output_densities", "Fukui_plus.cube"), comment="Fukui f+")

    # Compute the f- Fukui function
    f_minus = cubes[charge] - cubes[charge-1]
    f_minus.save(join("./output_densities", "Fukui_minus.cube"), comment="Fukui f-")

    # Compute the f0 Fukui function
    f_zero = cubes[charge+1] - cubes[charge-1]
    f_zero = f_zero.scale(0.5)
    f_zero.save(join("./output_densities", "Fukui_zero.cube"), comment="Fukui f0")
