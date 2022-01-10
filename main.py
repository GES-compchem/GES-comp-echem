import os
from jobdispatcher import JobDispatcher
from compechem import molecule
from compechem.calculators import tools, crest, xtb, orca

os.chdir("/mnt/f/Downloads")
xyz_folder = "xyz_files"


molecules = []

for file in os.listdir(xyz_folder):
    molecules.append(
        molecule.Molecule(os.path.join(xyz_folder, file), charge=0, spin=1)
    )

for mol in molecules:

    deprotomers = crest.deprotonate(mol, nproc=8, remove=True)

    for deprotomer in deprotomers:
        print(deprotomer.name, deprotomer.energies["total"])
