import os
from jobdispatcher import JobDispatcher
import compechem.calculators.calculators as cc

xyz_folder = "/mnt/f/Downloads/xyz_files"

molecules = []

for file in os.listdir(xyz_folder):
    molecules.append(cc.Molecule(os.path.join(xyz_folder, file), charge=0, spin=1))

# for molecule in molecules:
#     cc.deprotonate(molecule)
#     for deprotomer in molecule.deprotomers:
#         for line in deprotomer:
#             print(line)

for molecule in molecules:
    cc.opt_xtb(molecule, conformers=False, nproc=1)
    print(molecule.energies)
