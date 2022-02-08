import os, shutil
from tempfile import mkdtemp
from compechem.molecule import Molecule
from compechem.modules import tools


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
    print(f"INFO: {mol.name}, charge {mol.charge} spin {mol.spin} - CREST tautomer search")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_TAUT", dir=os.getcwd())

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
                Molecule(f"{mol.name}_tautomer_{num}.xyz", charge=mol.charge, spin=mol.spin,)
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
    print(f"INFO: {mol.name}, charge {mol.charge} spin {mol.spin} - CREST conformer search")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_CONF", dir=os.getcwd())

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
                Molecule(f"{mol.name}_conformer_{num}.xyz", charge=mol.charge, spin=mol.spin,)
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
    print(f"INFO: {mol.name}, charge {mol.charge} spin {mol.spin} - CREST deprotonation")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_DEPROT", dir=os.getcwd())

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
                Molecule(f"{mol.name}_deprotomer_{num}.xyz", charge=mol.charge - 1, spin=mol.spin,)
            )
            num += 1

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return deprotomers


def qcg_grow(
    solute, solvent, charge=None, spin=None, method="gfnff", nsolv=None, nproc=1, remove_tdir=True
):
    """Quantum Cluster Growth using CREST.

    Parameters
    ----------
    solute : Molecule object
        solute molecule to use in the calculation
    solvent : Molecule object
        solvent molecule to use in the calculation
    charge : int, optional
            total charge of the molecule. Default is taken from the solute molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the solute molecule.
    method : str
        method for the geometry optimizations, by default gfnff
        Alternative options: gfn1, gfn2
    nsolv : int
        number of solvent molecules to add to the cluster, by default None.
        If a number is not specified, the program will keep adding solvent
        molecules until convergence is reached, or 150 molecules are added.
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True

    Returns
    -------
    cluster : Molecule object
        Molecule object containing the explicitly solvated input molecule
    """

    if charge is None:
        charge = solute.charge
    if spin is None:
        spin = solute.spin

    parent_dir = os.getcwd()
    print(f"INFO: {solute.name}, charge {charge} spin {spin} - CREST QCG")

    tdir = mkdtemp(prefix=solute.name + "_", suffix="_QCG", dir=os.getcwd())

    os.chdir(tdir)
    solute.write_xyz("solute.xyz")
    solvent.write_xyz("solvent.xyz")

    if nsolv is not None:
        os.system(
            f"crest solute.xyz --qcg solvent.xyz --nsolv {nsolv} --{method} --alpb water --charge {charge} --uhf {spin-1} --T {nproc} > output.out 2>> output.out"
        )

    else:
        os.system(
            f"crest solute.xyz --qcg solvent.xyz --{method} --alpb water --charge {charge} --uhf {spin-1} --T {nproc} > output.out 2>> output.out"
        )

    cluster = Molecule("grow/cluster.xyz", charge, spin)

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return cluster

