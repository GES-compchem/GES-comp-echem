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

    mol1 = Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)

    xtb.opt(mol1)
    tools.info(mol1)

    mol2 = crest.deprotonate(mol1, nproc=ncores)[0]
    xtb.opt(mol2)
    tools.info(mol2)
    mol2.write_xyz("deprot.xyz")

    pka = algorithms.calculate_pka(mol1, mol2, method="gfn2")
    print(f"pKa: {pka}")
