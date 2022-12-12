import os, copy, shutil, sh
from typing import Dict, Tuple
from tempfile import mkdtemp
from compechem.config import get_ncores
from compechem.systems import Ensemble, System
from compechem.tools import process_output, process_density

import logging

logger = logging.getLogger(__name__)


class OrcaInput:
    """Interface for running Orca calculations."""

    def __init__(
        self,
        method: str = "PBE",
        basis_set: str = "def2-TZVP",
        aux_basis: str = "def2/J",
        solvent: str = None,
        optionals: str = "",
        MPI_FLAGS: str = "",
        ORCADIR: str = "$ORCADIR",
    ) -> None:
        """
        Parameters
        ----------
        method : str
            level of theory, by default "PBE"
        basis_set : str, optional
            basis set, by default "def2-TZVP"
        aux_basis : str, optional
            auxiliary basis set for RIJCOSX, by default "def2/J"
        solvent : str, optional
            SMD solvent, by default None
        optionals : str, optional
            optional keywords, by default ""
        MPI_FLAGS: str, optional
            string containing optional flags to be passed to MPI when launching an ora job.
            (e.g. `--bind-to-none` or `--use-hwthread-cpus`), by default ""
        ORCADIR: str, optional
            the path or environment variable containing the path to the ORCA folder, by
            default "$ORCADIR"
        """

        self.method = method
        self.basis_set = basis_set
        self.aux_basis = aux_basis
        self.solvent = solvent
        self.optionals = optionals
        self.__MPI_FLAGS = MPI_FLAGS
        self.__ORCADIR = ORCADIR

    def write_input(
        self,
        mol: System,
        job_info: Dict,
    ) -> None:

        mol.write_xyz(f"{mol.name}.xyz")

        input = (
            "%pal\n"
            f"  nprocs {job_info['ncores']}\n"
            "end\n\n"
            f"%maxcore {job_info['maxcore']}\n\n"
            f"! {self.method} {self.basis_set} {self.optionals}\n"
        )
        if self.aux_basis:
            input += f"! RIJCOSX {self.aux_basis}\n"

        if job_info["type"] == "spe":
            input += "\n"

        elif job_info["type"] == "opt":
            if self.solvent:
                input += "! Opt NumFreq\n\n"
            else:
                input += "! Opt Freq\n\n"

        elif job_info["type"] == "freq":
            if self.solvent:
                input += "! NumFreq\n\n"
            else:
                input += "! Freq\n\n"

        elif job_info["type"] == "nfreq":
            input += "! NumFreq\n\n"

        elif job_info["type"] == "scan":
            input += "! Opt\n"
            input += "%geom\n" "  scan\n" f"    {job_info['scan']}\n" "  end\n"
            if job_info["constraints"]:
                input += (
                    "  constraints\n" f"    {{ {job_info['constraints']} C }}\n" "  end\n"
                )
            if job_info["invertconstraints"]:
                input += "  invertConstraints true\n"
            input += "end\n\n"

        if self.solvent:
            input += "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n\n"

        if job_info["save_cubes"]:
            input += "%plots\n"
            input += "  Format Gaussian_Cube\n"
            input += f"  dim1 {job_info['cube_dim']}\n"
            input += f"  dim2 {job_info['cube_dim']}\n"
            input += f"  dim3 {job_info['cube_dim']}\n"
            input += '  ElDens("eldens.cube");\n'
            if job_info["spin"] != 1:
                input += '  SpinDens("spindens.cube");\n'
            input += "end\n\n"

        input += f"* xyzfile {job_info['charge']} {job_info['spin']} {mol.name}.xyz\n"

        with open("input.inp", "w") as inp:
            inp.writelines(input)

        return

    def spe(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
        save_cubes: bool = False,
        cube_dim: int = 250,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Single point energy calculation.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 350
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        save_cubes: bool, optional
            if set to True, will save a cube file containing electronic and spin densities,
            by default False.
        cube_dim: int, optional
            resolution for the cube files (default 250)
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output molecule containing the new energies.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} SPE")
        logger.debug(f"Running ORCA calculation on {ncores} cores and {maxcore} MB of RAM")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_spe",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "spe",
                    "ncores": ncores,
                    "maxcore": maxcore,
                    "charge": charge,
                    "spin": spin,
                    "save_cubes": save_cubes,
                    "cube_dim": cube_dim,
                },
            )

            os.system(f"{self.__ORCADIR}/orca input.inp > output.out {self.__MPI_FLAGS}")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])

            mulliken_charges, mulliken_spin_populations = parse_mulliken("output.out")

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

                newmol.properties = copy.copy(mol.properties)

                newmol.properties.add(self.method)
                newmol.properties[self.method].electronic_energy = electronic_energy

                newmol.properties[self.method].mulliken_charges = mulliken_charges
                newmol.properties[
                    self.method
                ].mulliken_spin_populations = mulliken_spin_populations

            else:
                mol.properties.add(self.method)
                mol.properties[self.method].electronic_energy = electronic_energy

                mol.properties[self.method].mulliken_charges = mulliken_charges
                mol.properties[
                    self.method
                ].mulliken_spin_populations = mulliken_spin_populations

            process_output(mol, self.method, "spe", charge, spin)

            if save_cubes:
                process_density(mol, self.method, "spe", charge, spin)

            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def opt(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
        save_cubes: bool = False,
        cube_dim: int = 250,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Geometry optimization + frequency analysis.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 350
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        save_cubes: bool, optional
            if set to True, will save a cube file containing electronic and spin densities,
            by default False.
        cube_dim: int, optional
            resolution for the cube files (default 250)
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output molecule containing the new geometry and energies.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} OPT")
        logger.debug(f"Running ORCA calculation on {ncores} cores and {maxcore} MB of RAM")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_opt",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "opt",
                    "ncores": ncores,
                    "maxcore": maxcore,
                    "charge": charge,
                    "spin": spin,
                    "save_cubes": save_cubes,
                    "cube_dim": cube_dim,
                },
            )

            os.system(f"{self.__ORCADIR}/orca input.inp > output.out {self.__MPI_FLAGS}")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])
                    if "G-E(el)" in line:
                        vibronic_energy = float(line.split()[-4])

            mulliken_charges, mulliken_spin_populations = parse_mulliken("output.out")

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

                newmol.properties.add(self.method)
                newmol.properties[self.method].electronic_energy = electronic_energy
                newmol.properties[self.method].vibronic_energy = vibronic_energy

                newmol.properties[self.method].mulliken_charges = mulliken_charges
                newmol.properties[
                    self.method
                ].mulliken_spin_populations = mulliken_spin_populations

            else:
                mol.load_xyz(f"{mol.name}.xyz")

                mol.properties.add(self.method)
                mol.properties[self.method].electronic_energy = electronic_energy
                mol.properties[self.method].vibronic_energy = vibronic_energy

                mol.properties[self.method].mulliken_charges = mulliken_charges
                mol.properties[
                    self.method
                ].mulliken_spin_populations = mulliken_spin_populations

            process_output(mol, self.method, "opt", charge, spin)

            if save_cubes:
                process_density(mol, self.method, "opt", charge, spin)

            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def freq(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Frequency analysis (analytical frequencies).

        Note: if the SMD solvation model is detected, defaults to numerical frequencies
        (analytical frequencies are not currently supported)

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 350
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output molecule containing the new energies.
        """

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} FREQ")
        logger.debug(f"Running ORCA calculation on {ncores} cores and {maxcore} MB of RAM")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_freq",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "freq",
                    "ncores": ncores,
                    "maxcore": maxcore,
                    "charge": charge,
                    "spin": spin,
                },
            )

            os.system(f"{self.__ORCADIR}/orca input.inp > output.out {self.__MPI_FLAGS}")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])
                    if "G-E(el)" in line:
                        vibronic_energy = float(line.split()[-4])

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

                newmol.properties = copy.copy(mol.properties)

                newmol.properties.add(self.method)
                newmol.properties[self.method].electronic_energy = electronic_energy
                newmol.properties[self.method].vibronic_energy = vibronic_energy

            else:

                mol.properties.add(self.method)
                mol.properties[self.method].electronic_energy = electronic_energy
                mol.properties[self.method].vibronic_energy = vibronic_energy

            process_output(mol, self.method, "freq", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def nfreq(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Frequency analysis (numerical frequencies).

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 350
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output molecule containing the new energies.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} NFREQ")
        logger.debug(f"Running ORCA calculation on {ncores} cores and {maxcore} MB of RAM")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_nfreq",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "nfreq",
                    "ncores": ncores,
                    "maxcore": maxcore,
                    "charge": charge,
                    "spin": spin,
                },
            )

            os.system(f"{self.__ORCADIR}/orca input.inp > output.out {self.__MPI_FLAGS}")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])
                    if "G-E(el)" in line:
                        vibronic_energy = float(line.split()[-4])

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

                newmol.properties = copy.copy(mol.properties)

                newmol.properties.add(self.method)
                newmol.properties[self.method].electronic_energy = electronic_energy
                newmol.properties[self.method].vibronic_energy = vibronic_energy

            else:

                mol.properties.add(self.method)
                mol.properties[self.method].electronic_energy = electronic_energy
                mol.properties[self.method].vibronic_energy = vibronic_energy

            process_output(mol, self.method, "numfreq", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def scan(
        self,
        mol: System,
        scan: str = None,
        constraints: str = None,
        invertconstraints: bool = False,
        ncores: int = None,
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
        remove_tdir: bool = True,
    ):
        """Relaxed surface scan.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        scan : str
            string for the scan section in the %geom block
        constraints : str
            string for the constraints section in the %geom block
        invertconstraints : bool, optional
            if True, treats the constraints block as the only coordinate NOT to constrain
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 350
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        scan_list : Ensemble object
            Output Ensemble containing the scan frames.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} SCAN")
        logger.debug(f"Running ORCA calculation on {ncores} cores and {maxcore} MB of RAM")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_scan",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "scan",
                    "ncores": ncores,
                    "maxcore": maxcore,
                    "charge": charge,
                    "spin": spin,
                    "scan": scan,
                    "constraints": constraints,
                    "invertconstraints": invertconstraints,
                },
            )

            os.system(f"{self.__ORCADIR}/orca input.inp > output.out {self.__MPI_FLAGS}")

            xyz_list = [
                xyz
                for xyz in os.listdir(".")
                if os.path.splitext(xyz)[1] == ".xyz"
                and os.path.splitext(xyz)[0][:5] == "input"
                and xyz != "input.xyz"
                and xyz != "input_trj.xyz"
            ]

            mol_list = []

            for xyz in xyz_list:
                index = xyz.split(".")[1]
                shutil.move(f"input.{index}.xyz", f"{mol.name}.{index}.xyz")
                mol_list.append(System(f"{mol.name}.{index}.xyz", charge=charge, spin=spin))

            ensemble = Ensemble(mol_list)

            process_output(mol, self.method, "scan", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            return ensemble


class M06(OrcaInput):
    def __init__(self):
        super().__init__(
            method="M062X",
            basis_set="def2-TZVP",
            aux_basis="def2/J",
            solvent="water",
            optionals="",
        )


class r2SCAN(OrcaInput):
    def __init__(self):
        super().__init__(
            method="r2SCAN-3c",
            basis_set="",
            aux_basis=None,
            solvent="water",
            optionals="",
        )


class CCSD(OrcaInput):
    def __init__(self):
        super().__init__(
            method="DLPNO-CCSD",
            basis_set="Extrapolate(2/3,ANO)",
            aux_basis="AutoAux",
            solvent="water",
            optionals="",
        )


def parse_mulliken(output: str) -> Tuple[Dict[int, float], Dict[int, float]]:
    """
    Given an ORCA output file the function will extract a dictionary containing the mulliken
    charges associated to each atom.

    Parameters
    ----------
    output: str
        The path to the ORCA output file.

    Raises
    ------
    ValueError
        Exception raised if the given path to the output file is not valid.

    Returns
    -------
    Dict[int, float]
        The dictionary containing the index of the atom as an integer key associated to the
        correspondent floating point Mulliken atomic charge value.
    Dict[int, float]
        The dictionary containing the index of the atom as an integer key associated to the
        correspondent floating point Mulliken atomic spin population value.
    """

    if os.path.exists(output):

        counter = 0
        charges, spins = {}, {}
        spin_available = False
        with open(output, "r") as file:

            # Count the number of "MULLIKEN ATOMIC CHARGES" sections in the file
            sections = file.read().count("MULLIKEN ATOMIC CHARGES")

            # Trace back to the beginning of the file
            file.seek(0)

            # Cycle over all the lines of the fuke
            for line in file:

                # If a "MULLIKEN ATOMIC CHARGES" section is found, increment the counter
                if "MULLIKEN ATOMIC CHARGES" in line:
                    counter += 1

                # If the index of the "MULLIKEN ATOMIC CHARGES" correspond with the last one
                # proceed with the file parsing else continue
                if counter == sections:

                    # Check if the section contains also the "SPIN POPULATIONS" column
                    if "SPIN POPULATIONS" in line:
                        spin_available = True

                    _ = file.readline()  # Skip the table line

                    # Iterate over the whole section reading line by line
                    while True:
                        buffer = file.readline()
                        if "Sum of atomic charges" in buffer:
                            break
                        else:
                            data = buffer.replace(":", "").split()
                            charges[int(data[0])] = float(data[2])

                            if spin_available:
                                spins[int(data[0])] = float(data[3])
                            else:
                                spins[int(data[0])] = 0
                else:
                    continue

                # If break has been called after mulliken has been modified the section end
                # has been reached, as such, break also from the reading operation
                if charges != {}:
                    break

    else:
        raise ValueError(f"File {output} not found!")

    return charges, spins
