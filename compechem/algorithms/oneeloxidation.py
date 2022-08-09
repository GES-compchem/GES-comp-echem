import os, pickle
import numpy as np

from compechem.config import get_ncores
from compechem.systems import System
from compechem.wrappers.xtb import XtbInput
from compechem.wrappers import crest
from compechem.functions.reorderenergies import reorder_energies
from compechem.functions.pka import calculate_pka
from compechem.functions.potential import calculate_potential

import logging

from typing import Iterator, Any

logger = logging.getLogger(__name__)

xtb = XtbInput()


class Species:
    """Class containing the singlets and radicals lists"""

    def __init__(self):
        self.singlets: list = []
        self.radicals: list = []


def find_highest_protonation_state(
    mol: System,
    method,
    ncores: int = None,
    maxcore: int = 350,
    conformer_search: bool = True,
):
    """Calculates the highest protonation state for a given molecule, as the first protomer
    with pKa < 0

    Parameters
    ----------
    mol : System object
        System object of the molecule under examination
    method : calculator
        Calculator object (i.e., XtbInput/OrcaInput object) giving the level of theory at which
        to evaluate the pKa for the protomers.
    ncores : int
        number of cores, by default all available cores
    maxcore : int, optional
        memory per core, in MB, by default 350
    conformer_search : Bool
        If True (default), also carries out conformer searches at all stages

    Returns
    -------
    currently_protonated : System object
        System object representing the highest protonation state for the input molecule
    """

    if ncores is None:
        ncores = get_ncores()

    if conformer_search:
        mol = crest.conformer_search(mol, ncores=ncores)[0]
    xtb.opt(mol, ncores=ncores, inplace=True)

    if type(method) != XtbInput:
        method.spe(mol, ncores=ncores, maxcore=maxcore, inplace=True)

    pKa = 0
    currently_deprotonated = mol

    while pKa > -2:

        protomer_list = crest.protonate(currently_deprotonated, ncores=ncores)

        if protomer_list:
            if type(method) != XtbInput:
                protomer_list = reorder_energies(
                    protomer_list,
                    ncores=ncores,
                    maxcore=maxcore,
                    method_opt=xtb,
                    method_el=method,
                    method_vib=xtb,
                )
            currently_protonated = protomer_list[0]

        # if protonation is unsuccessful (e.g., topology change), return the original molecule.
        else:
            return currently_deprotonated

        if conformer_search:
            currently_protonated = crest.conformer_search(
                currently_protonated, ncores=ncores
            )[0]

        xtb.opt(currently_protonated, inplace=True)
        if type(method) != XtbInput:
            method.spe(currently_protonated, ncores=ncores, maxcore=maxcore, inplace=True)

        try:
            pKa = calculate_pka(
                protonated=currently_protonated,
                deprotonated=currently_deprotonated,
                method_el=method.method,
                method_vib=xtb.method,
            )

        except:
            pKa = None

        logger.info(
            f"{currently_protonated.name} (charge {currently_protonated.charge} spin {currently_protonated.spin}) pKa = {pKa} ({method.method})"
        )

        if pKa is None:
            break
        else:
            currently_deprotonated = currently_protonated

    return currently_protonated


def calculate_deprotomers(
    mol: System,
    method,
    ncores: int = None,
    maxcore: int = 350,
    conformer_search: bool = True,
):
    """Calculates all the deprotomers with a pKa < 20 for a given System object

    Parameters
    ----------
    mol : System object
        System object of the molecule under examination
    method : calculator
        Calculator object (i.e., XtbInput/OrcaInput object) giving the level of theory at which
        to evaluate the pKa for the deprotomers.
    ncores : int
        number of cores, by default all available cores
    maxcore : int, optional
        memory per core, in MB, by default 350
    conformer_search : Bool
        If True (default), also carries out conformer searches at all stages

    Returns
    -------
    mol_list : list
        List containing all the deprotomers with pKa < 20 for the given molecule
    """

    if ncores is None:
        ncores = get_ncores()

    mol_list = []

    if conformer_search:
        mol = crest.conformer_search(mol, ncores=ncores)[0]
    xtb.opt(mol, ncores=ncores, inplace=True)

    if type(method) != XtbInput:
        method.spe(mol, ncores=ncores, maxcore=maxcore, inplace=True)

    pKa = 0
    i = 1
    currently_protonated = mol

    while pKa < 20:

        deprotomer_list = crest.deprotonate(currently_protonated, ncores=ncores)

        if deprotomer_list:
            if type(method) != XtbInput:
                deprotomer_list = reorder_energies(
                    deprotomer_list,
                    ncores=ncores,
                    maxcore=maxcore,
                    method_opt=xtb,
                    method_el=method,
                    method_vib=xtb,
                )
            currently_deprotonated = deprotomer_list[0]

        # if deprotonation is unsuccessful (e.g., topology change), save the molecule but with
        # sentinel values for pKa (which cannot be calculated)
        else:
            currently_protonated.properties.pka[method.method] = None
            mol_list.append(currently_protonated)
            break

        if conformer_search:
            currently_deprotonated = crest.conformer_search(
                currently_deprotonated, ncores=ncores
            )[0]

        xtb.opt(currently_deprotonated, inplace=True)
        if type(method) != XtbInput:
            method.spe(currently_deprotonated, ncores=ncores, maxcore=maxcore, inplace=True)

        try:
            pKa = calculate_pka(
                protonated=currently_protonated,
                deprotonated=currently_deprotonated,
                method_el=method.method,
                method_vib=xtb.method,
            )

        except:
            pKa = None

        logger.info(
            f"{currently_protonated.name} (charge {currently_protonated.charge} spin {currently_protonated.spin}) pKa{i} = {pKa} ({method.method})"
        )

        currently_protonated.properties.pka[method.method] = pKa

        mol_list.append(currently_protonated)

        if pKa is None:
            break
        else:
            currently_protonated = currently_deprotonated
            i += 1

    return mol_list


def generate_species(
    base_mol: System,
    method,
    ncores: int = None,
    maxcore: int = 350,
    conformer_search: bool = True,
    tautomer_search: bool = True,
):
    """Carries out all the calculations for the singlet and radical species to be used in the
    calculate_potential function, given a file path with a .xyz file
    
    Parameters
    ----------
    base_mol : System
        System object to be used in the calculation
    method : calculator
        Calculator object (i.e., XtbInput/OrcaInput object) giving the level of theory at which
        to evaluate the pKa for the deprotomers.
    ncores : int
        number of cores, by default all available cores
    maxcore : int, optional
        memory per core, in MB, by default 350
    conformer_search : Bool
        If True (default), also carries out a preliminary conformer search on the given structure
    tautomer_search : Bool
        If True (default), also carries out a preliminary conformer search on the given structure

    Returns
    -------
    species : Species
        Species object with the singlets and radicals for the given input molecule
    """

    if ncores is None:
        ncores = get_ncores()

    species = Species()

    molname = base_mol.name

    try:
        if conformer_search:
            base_mol = crest.conformer_search(
                base_mol, ncores=ncores, optionals="--noreftopo"
            )[0]
        if tautomer_search:
            tautomer_list = crest.tautomer_search(
                base_mol, ncores=ncores, optionals="--noreftopo"
            )
            if type(method) != XtbInput:
                tautomer_list = reorder_energies(
                    tautomer_list,
                    ncores=ncores,
                    maxcore=maxcore,
                    method_opt=xtb,
                    method_el=method,
                    method_vib=xtb,
                )
            base_mol = tautomer_list[0]

        ### SINGLET SPECIES ###
        singlet = xtb.spe(
            base_mol, ncores=ncores, charge=base_mol.charge, spin=base_mol.spin
        )
        singlet = find_highest_protonation_state(
            mol=singlet,
            method=method,
            ncores=ncores,
            maxcore=maxcore,
            conformer_search=conformer_search,
        )
        species.singlets = calculate_deprotomers(
            mol=singlet,
            method=method,
            ncores=ncores,
            maxcore=maxcore,
            conformer_search=conformer_search,
        )

        ### RADICAL SPECIES ###
        radical = xtb.spe(
            base_mol, ncores=ncores, charge=base_mol.charge + 1, spin=base_mol.spin + 1
        )
        radical = find_highest_protonation_state(
            mol=radical,
            method=method,
            ncores=ncores,
            maxcore=maxcore,
            conformer_search=conformer_search,
        )
        species.radicals = calculate_deprotomers(
            mol=radical,
            method=method,
            ncores=ncores,
            maxcore=maxcore,
            conformer_search=conformer_search,
        )

    except Exception as e:
        logger.error(f"Error occurred for {molname}! Skipping molecule")
        logger.error(e)
        return

    return species


def generate_potential_data(
    species: Species, method, pH_step: float = 1.0
) -> Iterator[Any]:
    """Calculates the 1-el oxidation potential for the given molecule in the pH range 0-14

    Parameters
    ----------
    species : Species
        Container object with the singlet and radical deprotomers for the input molecule
    method : calculator
        Calculator object (i.e., XtbInput/OrcaInput object) giving the level of theory at which
        to evaluate the pKa for the deprotomers.
    pH_step : float
        pH step at which the potential is calculated (by default, 1.0 pH units)


    Yields (generator)
    -------
    current_pH, potential

    """

    last_potential = None

    try:
        molname = species.singlets[0].name
        singlets = species.singlets
        radicals = species.radicals
    except AttributeError as e:
        logger.error(f"Missing suitable singlet/radical species in file {species}")
        logger.exception(e)

    logger.info(
        f"Generating potentials data for {species.singlets[0].name}, method: {method.method}, pH step: {pH_step}"
    )

    for current_pH in np.around(np.arange(0, 14 + pH_step, pH_step), 1):

        # getting deprotomers with pKa > current_pH or with sentinel value (deprotonation was
        # unsuccessful)
        for singlet in singlets:
            pka = singlet.properties.pka[method.method]
            if pka is None or pka > current_pH:
                current_singlet = singlet
                break

        for radical in radicals:
            pka = radical.properties.pka[method.method]
            if pka is None or pka > current_pH:
                current_radical = radical
                break

        try:
            potential = calculate_potential(
                oxidised=current_radical,
                reduced=current_singlet,
                pH=current_pH,
                method_el=method.method,
                method_vib=xtb.method,
            )
        except KeyError as e:
            logger.error(
                f"Requested level of theory not available in provided data for {molname}"
            )
            logger.exception(e)

        if last_potential is not None and abs(potential - last_potential) > 0.3:
            logger.warning(
                f"{molname} Potential changed by {abs(potential - last_potential)} at pH {current_pH}"
            )

        last_potential = potential

        yield current_pH, potential


def one_electron_oxidation_potentials(
    molecule: System,
    method,
    ncores: int = None,
    maxcore: int = 350,
    conformer_search: bool = True,
    tautomer_search: bool = True,
    pH_step: float = 1.0,
):

    if ncores is None:
        ncores = get_ncores()

    logger.debug(
        f"Requested 1-el oxidation calculation on {ncores} cores with {maxcore} MB of RAM"
    )

    os.makedirs("pickle_files", exist_ok=True)

    basename = molecule.name

    species = generate_species(
        base_mol=molecule,
        method=method,
        ncores=ncores,
        maxcore=maxcore,
        conformer_search=conformer_search,
        tautomer_search=tautomer_search,
    )

    pickle.dump(species, open(f"pickle_files/{basename}.species", "wb"))

    try:
        data_generator = generate_potential_data(
            species=species, method=method, pH_step=pH_step
        )
    except Exception as e:
        logger.error(f"Could not calculate potential data for {basename}!")
        logger.exception(e)
        return []

    return data_generator
