import os, shutil
from tempfile import mkdtemp
from compechem import molecule
from compechem.calculators import tools


def tautomer_search(mol, nproc=1, remove_tdir=True):
    """Tautomer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True

    Returns
    -------
    tautomers : list
        list containing the found tautomers, in order of ascending energy
    """

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - CREST tautomer search")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_crestTAUT", dir=os.getcwd())

    os.chdir(tdir)
    mol.write_xyz("geom.xyz")

    os.system(
        f"crest geom.xyz --alpb water --charge {mol.charge} --uhf {mol.spin-1} --mquick --fstrict --tautomerize -T {nproc} > output.out 2>> output.out"
    )

    tools.cyclization_check(mol, "geom.xyz", "tautomers.xyz")

    molsize = mol.atomcount
    with open("tautomers.xyz", "r") as f:
        numlines = int(sum(1 for line in f))

    tautomers = []

    with open("tautomers.xyz", "r") as f:
        num = 1
        j = 0
        while j < numlines:
            i = 0
            with open(f"{mol.name}_tautomer_{num}.xyz", "w") as out:
                while i < molsize + 2:
                    out.write(f.readline())
                    i += 1
                    j += 1
            tautomers.append(
                molecule.Molecule(
                    f"{mol.name}_tautomer_{num}.xyz", charge=mol.charge, spin=mol.spin,
                )
            )
            num += 1

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return tautomers


def conformer_search(mol, nproc=1, remove_tdir=True):
    """Conformer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True

    Returns
    -------
    conformers : list
        list containing the found conformers, in order of ascending energy
    """

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - CREST conformer search")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_crestCONF", dir=os.getcwd())

    os.chdir(tdir)
    mol.write_xyz("geom.xyz")

    os.system(
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --mquick -T {nproc} > conformers.out 2>> conformers.out"
    )

    molsize = mol.atomcount
    with open("crest_conformers.xyz", "r") as f:
        numlines = int(sum(1 for line in f))

    conformers = []

    with open("crest_conformers.xyz", "r") as f:
        num = 1
        j = 0
        while j < numlines:
            i = 0
            with open(f"{mol.name}_conformer_{num}.xyz", "w") as out:
                while i < molsize + 2:
                    out.write(f.readline())
                    i += 1
                    j += 1
            conformers.append(
                molecule.Molecule(
                    f"{mol.name}_conformer_{num}.xyz", charge=mol.charge, spin=mol.spin,
                )
            )
            num += 1

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return conformers


def deprotonate(mol, nproc=1, remove_tdir=True):
    """Deprotomer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True

    Returns
    -------
    deprotomers : list
        list containing the found deprotomers, in order of ascending energy
    """

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name} - CREST deprotonation")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_crestDEPROT", dir=os.getcwd())

    os.chdir(tdir)
    mol.write_xyz("geom.xyz")

    os.system(
        f"crest geom.xyz --alpb water --charge {mol.charge} --uhf {mol.spin-1} --deprotonate -T {nproc} > output.out 2>> output.out"
    )

    molsize = mol.atomcount - 1
    with open("deprotonated.xyz", "r") as f:
        numlines = int(sum(1 for line in f))

    j = 0

    deprotomers = []

    num = 1
    with open("deprotonated.xyz", "r") as f:
        while j < numlines:
            i = 0
            with open(f"{mol.name}_deprotomer_{num}.xyz", "w") as out:
                while i < molsize + 2:
                    out.write(f.readline())
                    i += 1
                    j += 1
            deprotomers.append(
                molecule.Molecule(
                    f"{mol.name}_deprotomer_{num}.xyz",
                    charge=mol.charge + 1,
                    spin=mol.spin,
                )
            )
            num += 1

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return deprotomers

