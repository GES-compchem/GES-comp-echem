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
        method="B97-D3",
        basis_set="def2-SVP",
        aux_basis="def2/J",
        nproc=ncores,
        maxcore=350,
        solvation=False,
        solvent="water",
        optionals="D3BJ",
    )

    xtb = XtbInput(
        method="gfn2",
        nproc=ncores,
        solvation=True,
        solvent="water",
        optionals="",
    )

    mol_gs = Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)
    xtb.opt(mol_gs)
    orca.spe(mol_gs)

    mol_1 = crest.deprotonate(mol_gs, nproc=ncores)[0]
    xtb.opt(mol_1)
    orca.spe(mol_1)

    mol_2 = crest.deprotonate(mol_1, nproc=ncores)[0]
    xtb.opt(mol_2)
    orca.spe(mol_2)

    tools.info(mol_gs, print_geometry=False)
    tools.info(mol_1, print_geometry=False)
    tools.info(mol_2, print_geometry=False)

    pka_1 = algorithms.calculate_pka(
        mol_gs, mol_1, method_el="B97-D3", method_vib="gfn2"
    )
    print(f"\npKa1: {pka_1}")

    pka_2 = algorithms.calculate_pka(
        mol_1, mol_2, method_el="B97-D3", method_vib="gfn2"
    )
    print(f"\npKa2: {pka_2}")

    # mol_rc = Molecule(os.path.join(xyz_folder, file), charge=1, spin=2)

    # xtb.opt(mol_gs)
    # xtb.opt(mol_rc)

    # mol_nr = crest.deprotonate(mol_rc, nproc=ncores)[0]
    # xtb.opt(mol_nr)

    # tools.info(mol_gs)
    # tools.info(mol_rc)
    # tools.info(mol_nr)

    # pka_rc = algorithms.calculate_pka(
    #     mol_rc, mol_nr, method_el="gfn2", method_vib="gfn2"
    # )
    # print(f"Radical cation pKa: {pka_rc}")

    # potential_nonPCET = algorithms.calculate_potential(
    #     oxidised=mol_rc,
    #     reduced=mol_gs,
    #     pH=3,
    #     method_el="gfn2",
    #     method_vib="gfn2",
    # )

    # print(f"Potential (non PCET): {potential_nonPCET} V")

    # potential_PCET = algorithms.calculate_potential(
    #     oxidised=mol_nr,
    #     reduced=mol_gs,
    #     pH=3,
    #     method_el="gfn2",
    #     method_vib="gfn2",
    # )

    # print(f"Potential (PCET): {potential_PCET} V")
