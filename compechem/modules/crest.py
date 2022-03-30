import os, copy
from tempfile import mkdtemp
from compechem.molecule import Molecule
from compechem.modules import tools


def tautomer_search(mol, nproc=1, remove_tdir=True, optionals=""):
    """Tautomer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True
    optionals : str, optional
        optional flags for calculation

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
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --mquick --fstrict --tautomerize {optionals} -T {nproc} > output.out 2>> output.err"
    )

    if os.path.exists("tautomers.xyz"):
        tautomers_to_check = tools.split_multixyz(mol, file="tautomers.xyz", suffix="t")

        tautomers = []

        while tautomers_to_check:
            tautomer = tautomers_to_check.pop(0)
            tautomer.write_xyz(f"{tautomer.name}.xyz")
            if tools.cyclization_check("geom.xyz", f"{tautomer.name}.xyz") is True:
                print(
                    f"WARNING: cyclization change spotted for {tautomer.name}. Removing from list."
                )
            else:
                tautomers.append(tautomer)

        tools.process_output(
            mol, "CREST", mol.charge, mol.spin, "tautomers", tdir, remove_tdir, parent_dir
        )
        return tautomers

    else:
        print(f"WARNING: no tautomers possible for {mol.name}. Ignoring tautomer search.")
        tools.process_output(
            mol, "CREST", mol.charge, mol.spin, "tautomers", tdir, remove_tdir, parent_dir
        )
        return [mol]


def conformer_search(mol, nproc=1, remove_tdir=True, optionals=""):
    """Conformer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True
    optionals : str, optional
        optional flags for calculation

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
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --mquick {optionals} -T {nproc} > output.out 2>> output.err"
    )

    if os.path.exists("crest_conformers.xyz"):
        conformers_to_check = tools.split_multixyz(mol, file="crest_conformers.xyz", suffix="c")

        conformers = []

        while conformers_to_check:
            conformer = conformers_to_check.pop(0)
            conformer.write_xyz(f"{conformer.name}.xyz")
            if tools.cyclization_check("geom.xyz", f"{conformer.name}.xyz") is True:
                print(
                    f"WARNING: cyclization change spotted for {conformer.name}. Removing from list."
                )
            else:
                conformers.append(conformer)

        tools.process_output(
            mol, "CREST", mol.charge, mol.spin, "conformers", tdir, remove_tdir, parent_dir
        )

        return conformers
    else:
        print("ERROR: conformer search failed.")
        os.chdir(parent_dir)
        return None


def deprotonate(mol, nproc=1, remove_tdir=True, optionals=""):
    """Deprotomer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True
    optionals : str, optional
        optional flags for calculation

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
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --deprotonate {optionals} -T {nproc} > output.out 2>> output.err"
    )

    if os.path.exists("deprotonated.xyz"):
        deprotomers_to_check = tools.split_multixyz(
            mol, file="deprotonated.xyz", suffix="d", charge=mol.charge - 1
        )

        deprotomers = []

        while deprotomers_to_check:
            deprotomer = deprotomers_to_check.pop(0)
            deprotomer.write_xyz(f"{deprotomer.name}.xyz")
            if tools.cyclization_check("geom.xyz", f"{deprotomer.name}.xyz") is True:
                print(
                    f"WARNING: cyclization change spotted for {deprotomer.name}. Removing from list."
                )
            else:
                deprotomers.append(deprotomer)

        tools.process_output(
            mol, "CREST", mol.charge, mol.spin, "deprotomers", tdir, remove_tdir, parent_dir
        )

        return conformers
    else:
        print("ERROR: deprotomer search failed.")
        os.chdir(parent_dir)
        return None

    return deprotomers


def protonate(mol, nproc=1, remove_tdir=True, optionals=""):
    """Protomer search using CREST.

    Parameters
    ----------
    mol : Molecule object
        input molecule to use in the calculation
    nproc : int, optional
        number of cores, by default 1
    remove_tdir : bool, optional
        temporary work directory will be removed, by default True
    optionals : str, optional
        optional flags for calculation

    Returns
    -------
    protomers : list
        list containing the found protomers, in order of ascending energy
    """

    parent_dir = os.getcwd()
    print(f"INFO: {mol.name}, charge {mol.charge} spin {mol.spin} - CREST protonation")

    tdir = mkdtemp(prefix=mol.name + "_", suffix="_PROT", dir=os.getcwd())

    os.chdir(tdir)
    mol.write_xyz("geom.xyz")

    os.system(
        f"crest geom.xyz --alpb water --chrg {mol.charge} --uhf {mol.spin-1} --protonate {optionals} -T {nproc} > output.out 2>> output.err"
    )

    if os.path.exists("protonated.xyz"):
        protomers_to_check = tools.split_multixyz(
            mol, file="protonated.xyz", suffix="p", charge=mol.charge + 1
        )

        protomers = []

        while protomers_to_check:
            protomer = protomers_to_check.pop(0)
            protomer.write_xyz(f"{protomer.name}.xyz")
            if tools.cyclization_check("geom.xyz", f"{protomer.name}.xyz") is True:
                print(
                    f"WARNING: cyclization change spotted for {protomer.name}. Removing from list."
                )
            else:
                protomers.append(deprotomer)

        tools.process_output(
            mol, "CREST", mol.charge, mol.spin, "protomers", tdir, remove_tdir, parent_dir
        )

        return conformers
    else:
        print("ERROR: protomer search failed.")
        os.chdir(parent_dir)
        return None

    return protomers


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
        os.chdir(parent_dir)
        return None

    tools.process_output(solute, "QCG", charge, spin, "grow", tdir, remove_tdir, parent_dir)

    return cluster


def qcg_ensemble(
    solute,
    solvent,
    charge=None,
    spin=None,
    method="gfn2",
    enslvl="gfn2",
    ensemble_choice="full_ensemble",
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
    ensemble_choice : str
        file containing the chosen ensemble after generation, by default "full_ensemble". Available
        options are:
            - "full_ensemble"
            - "final_ensemble"
            - "crest_best"
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
        ensemble = tools.split_multixyz(solute, file=f"ensemble/{ensemble_choice}.xyz", suffix="e")

    except:
        print("ERROR: cluster growth failed.")
        os.chdir(parent_dir)
        return None

    tools.process_output(solute, "QCG", charge, spin, "ensemble", tdir, remove_tdir, parent_dir)

    return ensemble
