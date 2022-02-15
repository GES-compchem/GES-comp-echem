import os, shutil
from tempfile import mkdtemp
from typing import final
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
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --mquick --fstrict --tautomerize -T {nproc} > output.out 2>> output.err"
    )

    tools.cyclization_check(mol, "geom.xyz", "tautomers.xyz")

    tautomers = tools.split_multixyz(mol, "tautomers.xyz")

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
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --mquick -T {nproc} > output.out 2>> output.err"
    )

    conformers = tools.split_multixyz(mol, "crest_conformers.xyz")

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
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --deprotonate -T {nproc} > output.out 2>> output.err"
    )

    deprotomers = tools.split_multixyz(mol, "deprotonated.xyz", charge=mol.charge - 1)

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return deprotomers


def qcg_grow(
    solute,
    solvent,
    charge=None,
    spin=None,
    method="gfn2",
    nsolv=0,
    nproc=1,
    optionals="",
    remove_tdir=True,
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
        method for the geometry optimizations, by default gfn2
        Alternative options: gfn1, gfnff
    nsolv : int
        number of solvent molecules to add to the cluster, by default 0 (unconstrained).
        If a number is not specified, the program will keep adding solvent
        molecules until convergence is reached, or 150 molecules are added.
    nproc : int, optional
        number of cores, by default 1
    optionals : str, optional
        optional flags for calculation
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
    print(
        f"INFO: {solute.name}, charge {charge} spin {spin} - CREST QCG GROW - {nsolv} solvent molecules"
    )

    tdir = mkdtemp(prefix=solute.name + "_", suffix="_QCG_G", dir=os.getcwd())

    os.chdir(tdir)
    solute.write_xyz("solute.xyz")
    solvent.write_xyz("solvent.xyz")

    os.system(
        f"crest solute.xyz --qcg solvent.xyz --nsolv {nsolv} --{method} --alpb water --chrg {charge} --uhf {spin-1} {optionals} --T {nproc} > output.out 2>> output.err"
    )

    solute.write_xyz(f"{solute.name}.xyz")
    cluster = Molecule(f"{solute.name}.xyz", charge, spin)

    try:
        cluster.update_geometry("grow/cluster.xyz")
    except:
        print("ERROR: cluster growth failed.")
        return

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return cluster


def qcg_ensemble(
    solute,
    solvent,
    charge=None,
    spin=None,
    method="gfn2",
    enslvl="gfn2",
    nsolv=0,
    nproc=1,
    optionals="",
    remove_tdir=True,
):
    """Quantum Cluster Growth + ensemble generation using CREST.

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
        method for the geometry optimizations, by default gfn2
        Alternative options: gfn1, gfnff
    enslvl : str
        method for the ensemble optimization, by default gfn2
        Alternative options: gfn1, gfnff
    nsolv : int
        number of solvent molecules to add to the cluster, by default 0 (unconstrained).
        If a number is not specified, the program will keep adding solvent
        molecules until convergence is reached, or 150 molecules are added.
    nproc : int, optional
        number of cores, by default 1
    optionals : str, optional
        optional flags for calculation
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True

    Returns
    -------
    cluster : Molecule object
        Molecule object containing the explicitly solvated input molecule, with updated energy
        coming from enseble generation (electronic contribution only). The vibronic contribution
        is taken from the input solute molecule (if present), while the electronic contribution
        is taken as the weighted average of all generated ensembles.
    """

    if charge is None:
        charge = solute.charge
    if spin is None:
        spin = solute.spin

    parent_dir = os.getcwd()
    print(
        f"INFO: {solute.name}, charge {charge} spin {spin} - CREST QCG ENSEMBLE - {nsolv} solvent molecules"
    )

    tdir = mkdtemp(prefix=solute.name + "_", suffix="_QCG_E", dir=os.getcwd())

    os.chdir(tdir)
    solute.write_xyz("solute.xyz")
    solvent.write_xyz("solvent.xyz")

    os.system(
        f"crest solute.xyz --qcg solvent.xyz --nsolv {nsolv} --{method} --ensemble --enslvl {enslvl} --alpb water --chrg {charge} --uhf {spin-1} {optionals} --T {nproc} > output.out 2>> output.err"
    )

    try:
        ensemble = tools.split_multixyz(solute, "ensemble/final_ensemble.xyz")

    except:
        print("ERROR: cluster growth failed.")
        os.chdir(parent_dir)
        return

    if remove_tdir is True:
        shutil.rmtree(tdir)
    os.chdir(parent_dir)

    return ensemble
