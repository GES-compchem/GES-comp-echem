import os, shutil
from tempfile import mkdtemp
from compechem.calculators import tools


def spe_ccsd(mol, nproc=1, maxcore=7500, remove=True):

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - CCSD SPE")

    wdir = mkdtemp(prefix=mol.name + "_", suffix="_ccsdSPE", dir=os.getcwd())

    os.chdir(wdir)
    tools.write_xyz(mol, "geom.xyz")

    with open("input.inp", "w") as inp:
        inp.write(
            f"%pal nproc {nproc} end\n"
            f"%maxcore {maxcore}\n"
            "! DLPNO-CCSD ano-pVTZ\n"
            "! RIJCOSX AutoAux\n"
            "%CPCM\n"
            "  SMD True\n"
            '  SMDsolvent "water"\n'
            "end\n"
            f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n"
        )

    os.system("$ORCADIR/orca input.inp > output.out")

    with open("output.out", "r") as out:
        for line in out:
            if "FINAL SINGLE POINT ENERGY" in line:
                mol.energies["electronic"] = float(line.split()[-1])
                mol.energies["total"] = (
                    mol.energies["electronic"] + mol.energies["vibronic"]
                )

    if remove is True:
        shutil.rmtree(wdir)
    os.chdir(parent_dir)


def spe_b97(mol, nproc=1, maxcore=350, remove=True):

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - B97 SPE")

    wdir = mkdtemp(prefix=mol.name + "_", suffix="_b97SPE", dir=os.getcwd())

    os.chdir(wdir)
    tools.write_xyz(mol, "geom.xyz")

    with open("input.inp", "w") as inp:
        inp.write(
            f"%pal nproc {nproc} end\n"
            f"%maxcore {maxcore}\n"
            "! B97-D3 D3BJ def2-TZVP\n"
            "! RIJCOSX def2/J\n"
            "%CPCM\n"
            "  SMD True\n"
            '  SMDsolvent "water"\n'
            "end\n"
            f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n"
        )

    os.system("$ORCADIR/orca input.inp > output.out")

    with open("output.out", "r") as out:
        for line in out:
            if "FINAL SINGLE POINT ENERGY" in line:
                mol.energies["electronic"] = float(line.split()[-1])
                mol.energies["total"] = (
                    mol.energies["electronic"] + mol.energies["vibronic"]
                )

    if remove is True:
        shutil.rmtree(wdir)
    os.chdir(parent_dir)


def opt_m06(mol, nproc=1, solvent=False, maxcore=350, remove=True):

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - M06-2X-D3 OPT")

    wdir = mkdtemp(prefix=mol.name + "_", suffix="_m06OPT", dir=os.getcwd())

    os.chdir(wdir)
    tools.write_xyz(mol, "geom.xyz")

    with open("input.inp", "w") as inp:
        inp.write(
            f"%pal nproc {nproc} end\n"
            f"%maxcore {maxcore}\n"
            "! M062X D3ZERO def2-TZVP Opt\n"
            "! RIJCOSX def2/J\n"
            "%geom\n"
            "  maxiter 500\n"
            "end\n"
        )
        if solvent is True:
            inp.write("%CPCM\n" "  SMD True\n" '  SMDsolvent "water"\n' "end\n")
        inp.write(f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n")

    os.system("$ORCADIR/orca input.inp > output.out")

    mol.update_geometry

    if remove is True:
        shutil.rmtree(wdir)
    os.chdir(parent_dir)


def freq_m06(mol, nproc=1, solvent=True, maxcore=350, remove=True):

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - M06-2X-D3 FREQ")

    wdir = mkdtemp(prefix=mol.name + "_", suffix="_m06FREQ", dir=os.getcwd())

    os.chdir(wdir)
    tools.write_xyz(mol, "geom.xyz")

    with open("input.inp", "w") as inp:
        inp.write(
            f"%pal nproc {nproc} end\n"
            f"%maxcore {maxcore}\n"
            "! M062X D3ZERO def2-TZVP NumFreq\n"
            "! RIJCOSX def2/J\n"
        )
        if solvent is True:
            inp.write("%CPCM\n" "  SMD True\n" '  SMDsolvent "water"\n' "end\n")
        inp.write(f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n")

    os.system("$ORCADIR/orca input.inp > output.out")

    with open("output.out", "r") as out:
        for line in out:
            if "FINAL SINGLE POINT ENERGY" in line:
                mol.energies["electronic"] = float(line.split()[-1])
                mol.energies["total"] = (
                    mol.energies["electronic"] + mol.energies["vibronic"]
                )

    if remove is True:
        shutil.rmtree(wdir)
    os.chdir(parent_dir)

