import os

# from jobdispatcher import JobDispatcher
from compechem.molecule import Molecule
from compechem.calculators import tools
from compechem.calculators.orca import OrcaInput
from compechem.calculators.xtb import XtbInput


# os.chdir("/mnt/f/Downloads")
os.chdir("/mnt/c/Users/LucaBabetto/Downloads")
xyz_folder = "xyz_files"

ncores = 4

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
        method="gfn2", nproc=ncores, solvation=True, solvent="water", optionals=""
    )

    mol_gs = Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)

    xtb.spe(mol_gs)
    orca.spe(mol_gs)
    tools.info(mol_gs)
    # print(mol_gs.energies["gfn2"].electronic)
    # print(mol_gs.energies["B97-D3"].electronic)

    # xtb.opt(mol_gs)
    # orca.spe(mol_gs)
    # print(mol_gs.energies)
    # print(mol_gs.energies["gfn2"].electronic)
    # print(mol_gs.energies["B97-D3"].electronic)
