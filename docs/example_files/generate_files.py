from shutil import rmtree
from compechem.systems import System
from compechem.engines.xtb import XtbInput

xtb = XtbInput()
water = System("./cursed_water.xyz")
water.name = "water"

xtb.opt(water, inplace=True)

water.save_json("water.json")

rmtree("error_files", ignore_errors=True)
rmtree("output_files", ignore_errors=True)
