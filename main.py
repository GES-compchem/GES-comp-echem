import os

# from jobdispatcher import JobDispatcher
from compechem.molecule import Molecule
from compechem.calculators import tools, crest, algorithms
from compechem.calculators.orca import OrcaInput
from compechem.calculators.xtb import XtbInput


# os.chdir("/mnt/f/Downloads")
os.chdir("/mnt/c/Users/LucaBabetto/Downloads")
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
    mol_rc = Molecule(os.path.join(xyz_folder, file), charge=1, spin=2)

    xtb.opt(mol_gs)
    xtb.opt(mol_rc)

    mol_nr = crest.deprotonate(mol_rc, nproc=ncores)[0]
    xtb.opt(mol_nr)

    tools.info(mol_gs)
    tools.info(mol_rc)
    tools.info(mol_nr)

    pka_rc = algorithms.calculate_pka(
        mol_rc, mol_nr, method_el="gfn2", method_vib="gfn2"
    )
    print(f"Radical cation pKa: {pka_rc}")

    potential_nonPCET = algorithms.calculate_potential(
        oxidised=mol_rc,
        reduced=mol_gs,
        pH=3,
        method_el="gfn2",
        method_vib="gfn2",
    )

    print(f"Potential (non PCET): {potential_nonPCET} V")

    potential_PCET = algorithms.calculate_potential(
        oxidised=mol_nr,
        reduced=mol_gs,
        pH=3,
        method_el="gfn2",
        method_vib="gfn2",
    )

    print(f"Potential (PCET): {potential_PCET} V")
