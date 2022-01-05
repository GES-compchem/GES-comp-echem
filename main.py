import os
from jobdispatcher import JobDispatcher
import compechem.calculators.calculators as cc

os.chdir("/mnt/f/Downloads")

xyz_folder = "./xyz_files"

molecules = []

for root, dirs, files in os.walk(os.path.abspath(xyz_folder)):
    for file in files:
        molecules.append(cc.Molecule(os.path.join(root, file), charge=0, spin=1))


for item in molecules:
    print(item.moleculename)

# pka_jobs = [calc(file, 0, 0) for file in xyz_files]

# jobprocessor = JobDispatcher([job.calculate_pka for job in pka_jobs],)

# results = jobprocessor.run()
