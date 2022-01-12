def calculate_pka(protonated, deprotonated, method_el, method_vib):
    """Calculates the pKa of a molecule, given the protonated and deprotonated
    forms.

    Parameters
    ----------
    protonated : Molecule object
        molecule in the protonated form
    deprotonated : Molecule object
        molecule in the deprotonated form
    method_el : str
        level of theory for the electronic energy
    method_vib : str
        level of theory for the vibronic contributions, by default None
        if method_vib is None, it defaults to the same level of theory of 
        method_el

    Returns
    -------
    pKa : float
        pKa of the molecule.
    """

    if method_vib is None:
        method_vib = method_el

    protonated_energy = (
        protonated.energies[method_el].electronic
        + protonated.energies[method_vib].vibronic
    )

    deprotonated_energy = (
        deprotonated.energies[method_el].electronic
        + deprotonated.energies[method_vib].vibronic
    )

    hydrogens = protonated.atomcount - deprotonated.atomcount

    if method_el == "gfn2":
        pka_correction = 164.22  # kcal/mol
    else:
        pka_correction = 0

    pka = -(
        (protonated_energy - deprotonated_energy) * 627.5
        + 270.29 * hydrogens
        - pka_correction * hydrogens
    ) / (2.303 * 1.98720425864083 / 1000 * 298.15)

    return pka

