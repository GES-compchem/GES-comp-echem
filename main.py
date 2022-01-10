import os
from jobdispatcher import JobDispatcher
from compechem import molecule
from compechem.calculators import tools, crest, xtb, orca

os.chdir("/mnt/f/Downloads")
xyz_folder = "xyz_files"


for file in os.listdir(xyz_folder):
    mol_gs = molecule.Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)
    xtb.opt(mol_gs)
    tools.info(mol_gs)

    mol_rc = molecule.Molecule(os.path.join(xyz_folder, file), charge=1, spin=2)
    xtb.opt(mol_rc)
    tools.info(mol_rc)

    mol_nr = crest.deprotonate(mol_rc)[0]
    xtb.opt(mol_nr)
    tools.info(mol_nr)

