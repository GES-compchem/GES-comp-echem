import os, shutil
from tempfile import mkdtemp
from compechem.calculators import tools


def opt(mol, nproc=1, remove=True):

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - xTB OPT")

    wdir = mkdtemp(prefix=mol.name + "_", suffix="_xtbOPT", dir=os.getcwd())

    os.chdir(wdir)
    tools.write_xyz(mol, "geom.xyz")

    os.system(
        f"xtb geom.xyz --alpb water --charge {mol.charge} --uhf {mol.spin-1} --ohess -P {nproc} > output.out 2>> output.out"
    )

    if tools.dissociation_check(mol) is True:
        return

    tools.cyclization_check(mol, "geom.xyz", "xtbopt.xyz")

    with open("output.out", "r") as out:
        for line in out:
            if "total free energy" in line:
                mol.energies["total"] = float(line.split()[-3])
            if "G(RRHO) contrib." in line:
                mol.energies["vibronic"] = float(line.split()[-3])

    mol.update_geometry("xtbopt.xyz")

    if remove is True:
        shutil.rmtree(wdir)
    os.chdir(parent_dir)


def spe(mol, nproc=1, remove=True):

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - xTB SPE")

    wdir = mkdtemp(prefix=mol.name + "_", suffix="_xtbSPE", dir=os.getcwd())

    os.chdir(wdir)
    tools.write_xyz(mol, "geom.xyz")

    os.system(
        f"xtb geom.xyz --alpb water --charge {mol.charge} --uhf {mol.spin-1} --hess -P {nproc} > output.out 2>> output.out"
    )

    with open("output.out", "r") as out:
        for line in out:
            if "total free energy" in line:
                mol.energies["total"] = float(line.split()[-3])
            if "G(RRHO) contrib." in line:
                mol.energies["vibronic"] = float(line.split()[-3])

    if remove is True:
        shutil.rmtree(wdir)
    os.chdir(parent_dir)
