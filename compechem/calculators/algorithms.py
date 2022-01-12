import os
from rdkit import Chem


def calculate_pka(protonated, deprotonated, method):
    """Calculates the pKa of a molecule, given the protonated and deprotonated
    forms.

    Parameters
    ----------
    protonated : Molecule object
        molecule in the protonated form
    deprotonated : Molecule object
        molecule in the deprotonated form
    method : str
        level of theory at which the pKa will be evaluated
        (necessary for applying corrections)

    Returns
    -------
    pKa : float
        pKa of the molecule.
    """

    protonated_energy = (
        protonated.energies[method].electronic
        + protonated.energies[method].vibronic
    )

    deprotonated_energy = (
        deprotonated.energies[method].electronic
        + deprotonated.energies[method].vibronic
    )

    print(f"protonated energy: {protonated_energy}")
    print(f"deprotonated energy: {deprotonated_energy}")

    hydrogens = protonated.atomcount - deprotonated.atomcount

    if method == "gfn2":
        pka_correction = 164.22  # kcal/mol
    else:
        pka_correction = 0

    pka = -(
        (protonated_energy - deprotonated_energy) * 627.5
        + 270.29 * hydrogens
        - pka_correction * hydrogens
    ) / (2.303 * 1.98720425864083 / 1000 * 298.15)

    return pka

