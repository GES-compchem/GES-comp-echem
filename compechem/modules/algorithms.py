from compechem.molecule import Molecule
from compechem.ensemble import Ensemble


def calculate_pka(protonated, deprotonated, method_el, method_vib=None):
    """Calculates the pKa of a molecule, given the protonated and deprotonated forms.

    Parameters
    ----------
    protonated : Molecule object
        molecule in the protonated form
    deprotonated : Molecule object
        molecule in the deprotonated form
    method_el : str
        level of theory for the electronic energy
    method_vib : str, optional
        level of theory for the vibronic contributions, if method_vib is None (default), 
        the pKa calculation will be carried out only considering the electronic component

    Returns
    -------
    pKa : float
        pKa of the molecule.
    """

    if type(protonated) or type(deprotonated) == Ensemble:
        print(f"ERROR: calculating pKa for Ensemble instead of Molecule. Currently not supported.")
        return None

    if protonated.atomcount - deprotonated.atomcount != 1:
        print(f"ERROR: {protonated.name} deprotomer differs for more than 1 atom.")
        return None

    if method_vib is None:
        protonated_energy = (protonated.energies[method_el].electronic) * 627.5  # kcal/mol
        deprotonated_energy = (deprotonated.energies[method_el].electronic) * 627.5  # kcal/mol

    else:
        protonated_energy = (
            protonated.energies[method_el].electronic + protonated.energies[method_vib].vibronic
        ) * 627.5  # kcal/mol

        deprotonated_energy = (
            deprotonated.energies[method_el].electronic + deprotonated.energies[method_vib].vibronic
        ) * 627.5  # kcal/mol

    proton_self_energy = 0

    if method_el == "gfn2":
        proton_self_energy = 164.22  # kcal/mol

    pka = ((deprotonated_energy + (-270.29 + proton_self_energy) - protonated_energy)) / (
        2.303 * 1.98720425864083 / 1000 * 298.15
    )

    return pka


def calculate_potential(oxidised, reduced, pH, method_el, method_vib=None):
    """Calculates the reduction potential of a molecule, given the oxidised and reduced forms.

    Parameters
    ----------
    oxidised : Molecule object
        molecule in the oxidised state
    reduced : Molecule object
        molecule in the reduced state
    pH : float
        pH at which the redox potential is evaluated
    method_el : str
        level of theory for the electronic energy
    method_vib : str, optional
        level of theory for the vibronic contributions, if method_vib is None (default), 
        the potential calculation will be carried out only considering the electronic component

    Returns
    -------
    potential : float
        reduction potential of the molecule.
    """

    if method_vib is None:
        oxidised_energy = (oxidised.energies[method_el].electronic) * 627.5  # kcal/mol
        reduced_energy = (reduced.energies[method_el].electronic) * 627.5  # kcal/mol

    else:
        oxidised_energy = (
            oxidised.energies[method_el].electronic + oxidised.energies[method_vib].vibronic
        ) * 627.5  # kcal/mol

        reduced_energy = (
            reduced.energies[method_el].electronic + reduced.energies[method_vib].vibronic
        ) * 627.5  # kcal/mol

    # note, count is reversed compared to pKa calculation
    hydrogens = reduced.atomcount - oxidised.atomcount

    reference_potential = 4.28  # V
    electron_self_energy = 0
    proton_self_energy = 0

    if method_el == "gfn2":
        electron_self_energy = 111.75  # kcal/mol
        proton_self_energy = 164.22  # kcal/mol

    potential = (
        -(
            reduced_energy
            + electron_self_energy
            - (oxidised_energy + (-270.29 + proton_self_energy) * hydrogens)
        )
        * (4184 / 96485)
        - reference_potential
        - 0.059 * pH * hydrogens
    )

    return potential

