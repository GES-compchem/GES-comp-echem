from compechem.molecule import System


def calculate_potential(
    oxidised: System,
    reduced: System,
    method_el: str,
    method_vib: str = None,
    pH: float = 7.0,
):
    """Calculates the reduction potential of a molecule, given the oxidised and reduced forms.

    Parameters
    ----------
    oxidised : System object
        molecule in the oxidised state
    reduced : System object
        molecule in the reduced state
    method_el : str
        level of theory for the electronic energy
    method_vib : str, optional
        level of theory for the vibronic contributions, if method_vib is None (default), 
        the potential calculation will be carried out only considering the electronic component
    pH : float
        pH at which the redox potential is evaluated, by default pH=7

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
