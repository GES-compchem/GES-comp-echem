import os, copy, shutil, sh
from tempfile import mkdtemp
from compechem.config import get_ncores
from compechem.systems import System, Energies
from compechem.systems import MDTrajectory
from compechem.tools import process_output
from compechem.tools import save_dftb_trajectory
from compechem.tools import compress_dftb_trajectory
from typing import Dict
import logging

logger = logging.getLogger(__name__)


class DFTBInput:
    """Interface for running DFTB+ calculations

    Attributes
    ----------
    method : str, optional
        level of theory, by default "DFTB". "xTB" also supported.
    parameters : str, optional
        parameters to be used for the DFTB Hamiltonian (by default 3ob)
    solver : str, optional
        LAPACK eigensolver method (check manual for available options)
    dispersion : bool, optional
        activates D3 dispersion corrections (off by default)
    parallel : str, optional
        selects either openmpi-parallel version (mpi) or shared memory version (nompi)
    verbose : bool, optional
        if set to True, saves the full DFTB+ output, otherwise, only the smaller files
    """

    def __init__(
        self,
        method: str = "DFTB",
        parameters: str = "3ob/3ob-3-1/",
        solver: str = None,
        thirdorder: bool = True,
        dispersion: bool = False,
        fermi: bool = False,
        fermi_temp: float = 300.0,
        parallel: str = "mpi",
        verbose: bool = True,
    ) -> None:
        """
        Parameters
        ----------
        method : str, optional
            level of theory, by default "DFTB". "xTB" also supported.
        parameters : str, optional
            parameters to be used for the DFTB Hamiltonian (by default 3ob)
        solver : str, optional
            LAPACK eigensolver method (check manual for available options)
        thirdorder : bool, optional
            activates the 3rd order terms in the DFTB Hamiltonian
        dispersion : bool, optional
            activates D3 dispersion corrections (off by default)
        fermi: bool, optional
            Fills the single particle levels according to a Fermi distribution (off by
            default).
        fermi_temp: float, optional
            Electronic temperature in Kelvin units. Note, this is ignored for thermostated
            simulations. By default, 300 K.
        parallel : str, optional
            selects either openmpi-parallel version (mpi) or shared memory version (nompi)
        verbose : bool, optional
            if set to True, saves the full DFTB+ output, otherwise, only the smaller files
        """

        self.method = method
        self.parameters = parameters
        self.solver = solver
        self.thirdorder = thirdorder
        self.dispersion = dispersion
        self.fermi = fermi
        self.fermi_temp = fermi_temp
        self.parallel = parallel
        self.verbose = verbose  # add to docs
        if self.verbose:
            self.output_path = "output.out"
        else:
            self.output_path = "/dev/null"

        self.atom_dict = {
            "Br": "d",
            "C": "p",
            "Ca": "p",
            "Cl": "d",
            "F": "p",
            "H": "s",
            "I": "d",
            "K": "p",
            "Mg": "p",
            "N": "p",
            "Na": "p",
            "O": "p",
            "P": "d",
            "S": "d",
            "Zn": "d",
        }

        self.hubbard_derivs = {
            "Br": -0.0573,
            "C": -0.1492,
            "Ca": -0.0340,
            "Cl": -0.0697,
            "F": -0.1623,
            "H": -0.1857,
            "I": -0.0433,
            "K": -0.0339,
            "Mg": -0.02,
            "N": -0.1535,
            "Na": -0.0454,
            "O": -0.1575,
            "P": -0.14,
            "S": -0.11,
            "Zn": -0.03,
        }

        self.spin_constants = {
            "H": [-0.072],
            "C": [-0.031, -0.025, -0.025, -0.023],
            "N": [-0.033, -0.027, -0.027, -0.026],
            "O": [-0.035, -0.030, -0.030, -0.028],
            "S": [-0.021, -0.017, 0.000, -0.017, -0.016, 0.000, 0.000, 0.000, -0.080],
        }

    def write_input(
        self,
        mol: System,
        job_info: Dict,
    ) -> None:

        mol.write_gen(f"{mol.name}.gen")

        with open(f"{mol.name}.gen") as file:
            lines = file.readlines()
            atom_types = lines[1].split()

        input = "Geometry = GenFormat {\n" f'  <<< "{mol.name}.gen"\n' "}\n\n"

        if job_info["type"] == "spe":
            input += "Driver = GeometryOptimization{\n" "  MaxSteps = 0\n" "}\n\n"

        elif job_info["type"] == "opt":
            input += (
                "Driver = GeometryOptimization{\n"
                f"  LatticeOpt = {'Yes' if job_info['latticeopt'] else 'No'}\n"
                "}\n\n"
            )

        elif job_info["type"] == "md_nvt":
            input += (
                "Driver = VelocityVerlet{\n"
                f"  TimeStep [fs] = {job_info['timestep']}\n"
                "  Thermostat = NoseHoover {\n"
                f"    Temperature [K] = {job_info['temperature']}\n"
                "    CouplingStrength [cm^-1] = 3200\n"
                "  }\n"
                f"  Steps = {job_info['steps']}\n"
                "  MovedAtoms = 1:-1\n"
                f"  MDRestartFrequency = {job_info['mdrestartfreq']}\n"
                "  Velocities [AA/ps] {\n"
            )
            for velocity in mol.velocities:
                input += f"    {velocity[1:]}"
            input += "  }\n" "}\n\n"

        elif job_info["type"] == "simulated_annealing":
            input += (
                "Driver = VelocityVerlet{\n"
                f"  TimeStep [fs] = {job_info['timestep']}\n"
                "  Thermostat = NoseHoover {\n"
                "    Temperature [Kelvin] = TemperatureProfile{\n"
                f"      constant 1 {job_info['start_temp']}\n"
                f"      linear {job_info['ramp_steps-1']} {job_info['target_temp']}\n"
                f"      constant {job_info['hold_steps']} {job_info['target_temp']}\n"
                f"      linear {job_info['ramp_steps']} {job_info['start_temp']}\n"
                "    }\n"
                "    CouplingStrength [cm^-1] = 3200\n"
                "  }\n"
                "  MovedAtoms = 1:-1\n"
                f"  MDRestartFrequency = {job_info['mdrestartfreq']}\n"
                "    Velocities [AA/ps] {\n"
            )
            for velocity in mol.velocities:
                input += f"    {velocity[1:]}"
            input += "  }\n" "}\n" "\n"

        input += (
            f"Hamiltonian = {self.method} {{\n"
            "  MaxSCCIterations = 500\n"
            f"  Charge = {job_info['charge']}\n"
        )

        if self.fermi:
            input += (
                "  Filling = Fermi {\n" f"    Temperature [K] = {self.fermi_temp}\n" "  }\n"
            )

        if job_info["spin"] != 1:
            input += (
                "  SpinPolarisation = Colinear {\n"
                f"    UnpairedElectrons = {job_info['spin']-1}\n"
                "  }\n"
                "  SpinConstants = {\n"
            )
            if self.method == "DFTB":
                input += "    ShellResolvedSpin = Yes\n"
            for atom in atom_types:
                input += (
                    f"    {atom} = {{\n"
                    f"      {' '.join(str(spin) for spin in self.spin_constants[atom])}\n"
                    "    }\n"
                )
            input += "  }\n"

        if self.method == "DFTB":
            if self.solver:
                input += f"  Solver = {self.solver} {{}}\n"
            input += (
                "  Scc = Yes\n"
                "  SlaterKosterFiles = Type2FileNames {\n"
                f'    Prefix = "{self.parameters}"\n'
                '    Separator = "-"\n'
                '    Suffix = ".skf"\n'
                "  }\n"
                "  MaxAngularMomentum {\n"
            )
            for atom in atom_types:
                input += f'    {atom} = "{self.atom_dict[atom]}"\n'
            input += "  }\n"
            if mol.periodic:
                input += "  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n"
            if self.thirdorder:
                input += "  ThirdOrderFull = Yes\n" "  HubbardDerivs {\n"
                for atom in atom_types:
                    input += f"    {atom} = {self.hubbard_derivs[atom]}\n"
                input += (
                    "  }\n" "  HCorrection = Damping {\n" "    Exponent = 4.00\n" "  }\n"
                )
            if self.dispersion:
                input += (
                    "  Dispersion = SimpleDftD3 {\n"
                    "    a1 = 0.746\n"
                    "    a2 = 4.191\n"
                    "    s6 = 1.0\n"
                    "    s8 = 3.209\n"
                    "  }\n"
                )
            input += "}\n"

        elif self.method == "xTB":
            if self.solver:
                input += f"  Solver = {self.solver} {{}}\n"
            self.parameters = "gfn2"
            input += '  Method = "GFN2-xTB"\n'
            if mol.periodic:
                input += "  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n"
            input += "}\n"

        input += "\n" "ParserOptions {\n" "  ParserVersion = 11\n" "}"

        with open("dftb_in.hsd", "w") as inp:
            inp.writelines(input)

        return

    def spe(
        self,
        mol: System,
        ncores: int = None,
        maxcore=None,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Single point energy calculation.

        Parameters
        ----------
        mol : System object
            Input molecule to use in the calculation.
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True

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
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

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
                    "charge": charge,
                    "spin": spin,
                },
            )

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > output.out 2>> output.err")

            elif self.parallel == "nompi":
                os.system(f"dftb+ > output.out 2>> output.err")

            with open("output.out", "r") as out:
                for line in out:
                    if "Total Energy" in line:
                        electronic_energy = float(line.split()[2])

            vibronic_energy = None

            if self.method in mol.energies:
                vibronic_energy = mol.energies[self.method].vibronic

            if inplace is False:

                mol.write_xyz(f"{mol.name}.xyz")

                newmol = System(f"{mol.name}.xyz", charge, spin, mol.box_side)

                newmol.energies = copy.copy(mol.energies)

                newmol.energies[self.method] = Energies(
                    method=self.method,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

            else:
                mol.energies[self.method] = Energies(
                    method=self.method,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

            process_output(mol, self.method, "spe", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def opt(
        self,
        mol: System,
        latticeopt: bool = False,
        ncores: int = None,
        maxcore=None,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Geometry optimization.

        Parameters
        ----------
        mol : System object
            Input molecule to use in the calculation.
        latticeopt : bool, optional
            If True, also optimize PBC conditions. By default, False
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True

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

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} OPT")
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

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
                    "charge": charge,
                    "spin": spin,
                    "latticeopt": latticeopt,
                },
            )

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > output.out 2>> output.err")

            elif self.parallel == "nompi":
                os.system(f"dftb+ > output.out 2>> output.err")

            with open("output.out", "r") as out:
                for line in out:
                    if "Total Energy" in line:
                        electronic_energy = float(line.split()[2])

            vibronic_energy = None

            if self.method in mol.energies:
                vibronic_energy = mol.energies[self.method].vibronic

            if inplace is False:

                mol.write_xyz(f"{mol.name}.xyz")
                newmol = System(f"{mol.name}.xyz", charge, spin, mol.box_side)

                newmol.energies = copy.copy(mol.energies)

                newmol.energies[self.method] = Energies(
                    method=self.method,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

                newmol.update_geometry("geo_end.xyz")

            else:
                mol.energies[self.method] = Energies(
                    method=self.method,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

                mol.update_geometry("geo_end.xyz")

            process_output(mol, self.method, "spe", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def md_nvt(
        self,
        mol: System,
        steps: int,
        timestep: float = 1.0,
        temperature: float = 298.0,
        mdrestartfreq: int = 100,
        box_side: float = None,
        ncores: int = None,
        maxcore=None,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
        compress_traj: bool = True,
    ):
        """Molecular Dynamics simulation in the Canonical Ensemble (NVT).

        Parameters
        ----------
        mol : System object
            Input molecule to use in the calculation.
        steps : int
            Total steps of the simulation
        timestep : float, optional
            Time step (in fs) for the simulation.
        temperature : float, optional
            Temperature (in Kelvin) of the simulation
        mdrestartfreq : int, optional
            MD information is printed to md.out every mdrestartfreq steps, by default 100
        box_side : float, optional
            for periodic systems, defines the length (in Å) of the box side
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True
        compress_traj : bool, optional
            if True, parses the geo.end and md.out files into a single, smaller file, which
            is then zipped in an archive.

        Returns
        -------
        trajectory : Ensemble object
            Ensemble containing the NVT MD trajectory data
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin
        if box_side is None:
            box_side = mol.box_side

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} NVT MD")
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_md_nvt",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "md_nvt",
                    "charge": charge,
                    "spin": spin,
                    "timestep": timestep,
                    "temperature": temperature,
                    "steps": steps,
                    "mdrestartfreq": mdrestartfreq,
                },
            )

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > {self.output_path} 2>> output.err")
            elif self.parallel == "nompi":
                os.environ["OMP_NUM_THREADS"] = f"{ncores}"
                os.system(f"dftb+ > {self.output_path} 2>> output.err")

            import random, string

            if inplace is False:
                newmol = System("geo_end.xyz", charge, spin, mol.box_side)
                newmol.energies = copy.copy(mol.energies)

            else:
                mol.update_geometry("geo_end.xyz")

            suffix = "".join(random.choices(string.ascii_letters + string.digits, k=4))

            if compress_traj:
                compress_dftb_trajectory(f"{mol.name}_{charge}_{spin}")
                os.makedirs("../MD_trajectories", exist_ok=True)
                shutil.move(
                    f"{mol.name}_{charge}_{spin}.zip",
                    f"../MD_trajectories/{mol.name}_{charge}_{spin}.zip",
                )

            save_dftb_trajectory(f"{mol.name}_{charge}_{spin}_{suffix}")

            if mol.periodic:
                with open(f"../MD_data/{mol.name}_{charge}_{spin}_{suffix}.pbc", "w") as f:
                    f.write(f"{mol.box_side}")

            process_output(mol, self.method, "md_nvt", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

        trajectory = MDTrajectory(f"{mol.name}_{charge}_{spin}_{suffix}", self.method)

        return trajectory

    def simulated_annealing(
        self,
        mol: System,
        start_temp: float = 1.0,
        target_temp: float = 2000.0,
        ramp_steps: int = 500,
        hold_steps: int = 1000,
        timestep: float = 1.0,
        mdrestartfreq: int = 100,
        box_side: float = None,
        ncores: int = None,
        maxcore=None,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
        compress_traj: bool = True,
    ):
        """Molecular Dynamics simulated annealing simulation in the Canonical Ensemble (NVT)

        Parameters
        ----------
        mol : System object
            Input molecule to use in the calculation.
        start_temp: float, optional
            Starting temperature (default, 1K)
        target_temp: float, optional
            Maximum temperature reached during the simulation (default, 2000K)
        ramp_steps: int, optional
            Number of MD steps for the heating/cooling ramps (default, 500 steps)
        hold_steps: int, optional
            Number of MD steps held at target_temp (default, 1000 steps)
        timestep : float, optional
            Time step (in fs) for the simulation.
        mdrestartfreq : int, optional
            MD information is printed to md.out every mdrestartfreq steps, by default 100
        box_side : float, optional
            for periodic systems, defines the length (in Å) of the box side
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True
        compress_traj : bool, optional
            if True, parses the geo.end and md.out files into a single, smaller file.

        Returns
        -------
        trajectory : Ensemble object
            Ensemble containing the NVT MD trajectory data
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin
        if box_side is None:
            box_side = mol.box_side

        logger.info(
            f"{mol.name}, charge {charge} spin {spin} - {self.method} Simulated Annealing"
        )
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")
        logger.debug(
            f"Heating/cooling between {start_temp}K and {target_temp}K for {ramp_steps} steps and holding max temp for {hold_steps} steps"
        )

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_anneal",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            self.write_input(
                mol=mol,
                job_info={
                    "type": "simulated_annealing",
                    "charge": charge,
                    "spin": spin,
                    "timestep": timestep,
                    "start_temp": start_temp,
                    "ramp_steps": ramp_steps,
                    "target_temp": target_temp,
                    "hold_steps": hold_steps,
                    "mdrestartfreq": mdrestartfreq,
                },
            )

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > {self.output_path} 2>> output.err")
            elif self.parallel == "nompi":
                os.environ["OMP_NUM_THREADS"] = f"{ncores}"
                os.system(f"dftb+ > {self.output_path} 2>> output.err")

            import random, string

            if inplace is False:
                newmol = System("geo_end.xyz", charge, spin, mol.box_side)
                newmol.energies = copy.copy(mol.energies)

            else:
                mol.update_geometry("geo_end.xyz")

            suffix = "".join(random.choices(string.ascii_letters + string.digits, k=4))

            if compress_traj:
                compress_dftb_trajectory(mol.name)
                os.makedirs("../MD_trajectories", exist_ok=True)
                shutil.move(f"{mol.name}.zip", f"../MD_trajectories/{mol.name}.zip")

            save_dftb_trajectory(f"{mol.name}_{suffix}")

            if mol.periodic:
                with open(f"../MD_data/{mol.name}_{suffix}.pbc", "w") as f:
                    f.write(f"{mol.box_side}")

            process_output(mol, self.method, "anneal", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

        trajectory = MDTrajectory(f"{mol.name}_{suffix}", self.method)

        return trajectory
