import os

# from jobdispatcher import JobDispatcher
from compechem.molecule import Molecule
from compechem.calculators import tools, crest, algorithms
from compechem.calculators.orca import OrcaInput
from compechem.calculators.xtb import XtbInput


os.chdir("/mnt/f/Downloads")
# os.chdir("/mnt/c/Users/LucaBabetto/Downloads")
xyz_folder = "xyz_files"

ncores = 8

for file in os.listdir(xyz_folder):

    orca = OrcaInput(
        method="M062X",
        basis_set="def2-TZVP",
        aux_basis="def2/J",
        nproc=ncores,
        maxcore=350,
        solvation=True,
        solvent="water",
        optionals="D3ZERO DEFGRID3",
    )

    xtb = XtbInput(
        method="gfn2",
        nproc=ncores,
        solvation=True,
        solvent="water",
        optionals="",
    )

    # CALCOLO PKA

    # mol_gs = Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)
    # xtb.opt(mol_gs)
    # orca.spe(mol_gs)

    # mol_1 = crest.deprotonate(mol_gs, nproc=ncores)[0]
    # xtb.opt(mol_1)
    # orca.spe(mol_1)

    # mol_2 = crest.deprotonate(mol_1, nproc=ncores)[0]
    # xtb.opt(mol_2)
    # orca.spe(mol_2)

    # tools.info(mol_gs, print_geometry=False)
    # tools.info(mol_1, print_geometry=False)
    # tools.info(mol_2, print_geometry=False)

    # pka_1 = algorithms.calculate_pka(
    #     mol_gs, mol_1, method_el="B97-D3", method_vib="gfn2"
    # )
    # print(f"\npKa1: {pka_1}")

    # pka_2 = algorithms.calculate_pka(
    #     mol_1, mol_2, method_el="B97-D3", method_vib="gfn2"
    # )
    # print(f"\npKa2: {pka_2}")

    # CALCOLO POTENZIALE

    mol_gs = Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)
    mol_rc = Molecule(os.path.join(xyz_folder, file), charge=1, spin=2)

    xtb.opt(mol_gs)
    orca.spe(mol_gs)
    xtb.opt(mol_rc)
    orca.spe(mol_rc)

    mol_nr = crest.deprotonate(mol_rc, nproc=ncores)[0]
    xtb.opt(mol_nr)
    orca.spe(mol_nr)

    tools.info(mol_gs)
    tools.info(mol_rc)
    tools.info(mol_nr)

    pka_rc1 = algorithms.calculate_pka(
        mol_rc, mol_nr, method_el="gfn2", method_vib="gfn2"
    )
    print(f"\nRadical cation pKa (gfn2): {pka_rc1}")

    pka_rc2 = algorithms.calculate_pka(
        mol_rc, mol_nr, method_el="M062X", method_vib="gfn2"
    )
    print(f"\nRadical cation pKa (M062X): {pka_rc2}")

    potential_nonPCET1 = algorithms.calculate_potential(
        oxidised=mol_rc,
        reduced=mol_gs,
        pH=1,
        method_el="gfn2",
        method_vib="gfn2",
    )

    print(f"\nPotential (non PCET) (gfn2): {potential_nonPCET1} V")

    potential_PCET1 = algorithms.calculate_potential(
        oxidised=mol_nr,
        reduced=mol_gs,
        pH=3,
        method_el="gfn2",
        method_vib="gfn2",
    )

    print(f"\nPotential (PCET) (gfn2): {potential_PCET1} V")

    potential_nonPCET2 = algorithms.calculate_potential(
        oxidised=mol_rc,
        reduced=mol_gs,
        pH=1,
        method_el="M062X",
        method_vib="gfn2",
    )

    print(f"\nPotential (non PCET) (M062X): {potential_nonPCET2} V")

    potential_PCET2 = algorithms.calculate_potential(
        oxidised=mol_nr,
        reduced=mol_gs,
        pH=3,
        method_el="M062X",
        method_vib="gfn2",
    )

    print(f"\nPotential (PCET) (M062X): {potential_PCET2} V")
