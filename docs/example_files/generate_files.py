# water opt
from compechem.systems import System
from compechem.engines.xtb import XtbInput

# xtb = XtbInput()
# water = System("./cursed_water.xyz")
# water.name = "water"

# xtb.opt(water, inplace=True)

# print("saving water.json")
# water.save_json("water.json")

# solvated cube
from compechem.wrappers.packmol import packmol_cube

solute = System("urea.xyz")
solvent = System("water.xyz")

solvated_cube = packmol_cube(solute, solvent, nsolv=50, target_dens=997)

print("saving solvated_cube.json")
solvated_cube.save_json("solvated_cube.json")

# generate acetaldehyde with Fukui functions
from compechem.engines.orca import OrcaInput
from compechem.functions.fukui import calculate_fukui

mol = System("./acetaldehyde.xyz")
orca = OrcaInput(method="PBE", basis_set="def2-SVP")

orca.opt(mol, inplace=True, ncores=4)
calculate_fukui(mol, orca, ncores=4)

print("saving acetaldehyde.json")
mol.save_json("acetaldehyde.json")

# remove tempfiles
from shutil import rmtree

rmtree("error_files", ignore_errors=True)
rmtree("output_files", ignore_errors=True)
rmtree("packmol_files", ignore_errors=True)
