from compechem.ensemble import Ensemble
from compechem.molecule import Molecule


def calculate_pka(
    protonated: Molecule, deprotonated: Molecule, method_el: str, method_vib: str = None
):
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

    if type(protonated) == Ensemble or type(deprotonated) == Ensemble:
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
