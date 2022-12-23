from os.path import join
from copy import deepcopy
from typing import List, Dict, Union

from compechem.systems import System
from compechem.engines.orca import OrcaInput
from compechem.tools.cubetools import Cube


FUKUI_CUBE_WARNING = "Fukui function cube file, the atomic charges column contains the localized Mulliken-charge-based fukui values"


def calculate_fukui(
    molecule: System,
    orca: OrcaInput,
    spins_states: Union[None, List[int]] = None,
    cube_dim: int = 250,
    ncores: int = None,
    maxcore: int = 1000,
) -> Dict[str, Dict[str, List[float]]]:
    """
    Computes the Fukui f+, f- and f0 functions starting from a given input molecule. The
    functions are saved in Gaussian cube compatible format and stored in the 
    `output_densities` folder. Localized Fukui functions are also computed from the Mulliken
    and Hirshfeld charges and saved in the molecule properties. Please notice that the charges
    in the `.fukui.cube` produced are replaced by the Mulliken condensed fukui functions.

    Parameters
    ----------
    molecule: System
        The System object containing the geometry of the selected molecule (if the geometry
        has not been optimized, please enable the optimize option)
    orca: OrcaInput
        The orca input wrapper object that defines the protocol to be used in the calculation.
    spin_states: Union[None, List[int]]
        If set to None, when adding or subtracting electrons will automatically switch the
        spin state from singlet to doublet and vice versa (Maximum one unpaired electrons).
        If manually set to `List[int]` will force the spin multeplicities according to the
        user specified values. The order of the spin multiplicity values is: molecule with
        one electron added (-1), the molecule as it is (0) and the molecule  with one of its
        electrons removed (+1).
    cube_dim: int
        Resolution for the cube files (default 250)
    ncores : int, optional
        The number of cores, by default all available cores (None)
    maxcore: int
        The maximum amount of memory in MB to be allocated for each core.
    """
    # RUN ALL THE REQUIRED CALCULATIONS
    # --------------------------------------------------------------------------------------
    # Compute a single point for the molecule
    orca.spe(
        molecule,
        save_cubes=True,
        cube_dim=cube_dim,
        hirshfeld=True,
        inplace=True,
        ncores=ncores,
        maxcore=maxcore,
    )

    # Compute a single point for the molecule with the addition of one electron.
    cation = deepcopy(molecule)
    cation.charge += 1

    if spins_states is not None:
        cation.spin = spins_states[2]
    else:
        cation.spin = 1 if molecule.spin == 2 else 2

    orca.spe(
        cation,
        save_cubes=True,
        cube_dim=cube_dim,
        hirshfeld=True,
        inplace=True,
        ncores=ncores,
        maxcore=maxcore
    )

    # Compute a single point for the molecule with the subtraction of one electron.
    anion = deepcopy(molecule)
    anion.charge -= 1

    if spins_states is not None:
        anion.spin = spins_states[0]
    else:
        anion.spin = 1 if molecule.spin == 2 else 2

    orca.spe(
        anion,
        save_cubes=True,
        cube_dim=cube_dim,
        hirshfeld=True,
        inplace=True,
        ncores=ncores,
        maxcore=maxcore
    )

    # COMPUTE CONDENSED FUKUI FUNCTIONS
    #---------------------------------------------------------------------------------------
    # Compute the localized fukui function values from the Mulliken charges
    localized_fukui_mulliken = {"f+": [], "f-": [], "f0": []}
    for atom in range(molecule.geometry.atomcount):
        localized_fukui_mulliken["f+"].append(
            -(anion.properties.mulliken_charges[atom] - molecule.properties.mulliken_charges[atom])
        )
        localized_fukui_mulliken["f-"].append(
            -(molecule.properties.mulliken_charges[atom] - cation.properties.mulliken_charges[atom])
        )
        localized_fukui_mulliken["f0"].append(
            -(anion.properties.mulliken_charges[atom] - cation.properties.mulliken_charges[atom]) / 2
        )
    
    molecule.properties.set_condensed_fukui_mulliken(localized_fukui_mulliken, orca)

    # Compute the localized fukui function values from the Hirshfeld charges
    localized_fukui_hirshfeld = {"f+": [], "f-": [], "f0": []}
    for atom in range(molecule.geometry.atomcount):
        localized_fukui_hirshfeld["f+"].append(
            -(anion.properties.hirshfeld_charges[atom] - molecule.properties.hirshfeld_charges[atom])
        )
        localized_fukui_hirshfeld["f-"].append(
            -(molecule.properties.hirshfeld_charges[atom] - cation.properties.hirshfeld_charges[atom])
        )
        localized_fukui_hirshfeld["f0"].append(
            -(anion.properties.hirshfeld_charges[atom] - cation.properties.hirshfeld_charges[atom]) / 2
        )
    
    molecule.properties.set_condensed_fukui_hirshfeld(localized_fukui_hirshfeld, orca)

    # LOAD THE VOLUMETRIC FILES AND COMPUTE THE FUKUI FUNCTIONS
    # --------------------------------------------------------------------------------------
    # Load cubes from the output_densities folder
    cubes: Dict[int, Cube] = {}
    for mol in [cation, molecule, anion]:
        cubename = f"{mol.name}_{mol.charge}_{mol.spin}_{orca.output_suffix}_spe.eldens.cube"
        cube = Cube.from_file(join("./output_densities", cubename))
        cubes[mol.charge] = cube

    # Check if all the densities have been loaded correctly
    if len(cubes) != 3:
        raise RuntimeError(f"Three cube files expected, {len(cubes)} found.")

    # Compute the f+ Fukui function
    f_plus = cubes[anion.charge] - cubes[molecule.charge]
    f_plus.charges = localized_fukui_mulliken["f+"]
    f_plus.save(
        join("./output_densities", f"{molecule.name}_Fukui_plus.fukui.cube"),
        comment_1st=FUKUI_CUBE_WARNING,
        comment_2nd="Fukui f+",
    )

    # Compute the f- Fukui function
    f_minus = cubes[molecule.charge] - cubes[cation.charge]
    f_minus.charges = localized_fukui_mulliken["f-"]
    f_minus.save(
        join("./output_densities", f"{molecule.name}_Fukui_minus.fukui.cube"),
        comment_1st=FUKUI_CUBE_WARNING,
        comment_2nd="Fukui f-",
    )

    # Compute the f0 Fukui function
    f_zero = (cubes[anion.charge] - cubes[cation.charge]).scale(0.5)
    f_zero.charges = localized_fukui_mulliken["f0"]
    f_zero.save(
        join("./output_densities", f"{molecule.name}_Fukui_zero.fukui.cube"),
        comment_1st=FUKUI_CUBE_WARNING,
        comment_2nd="Fukui f0",
    )
    