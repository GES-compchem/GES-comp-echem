import os, copy, shutil, sh
from tracemalloc import start
from tempfile import mkdtemp
from compechem.config import get_ncores
from compechem.systems import System, Energies
from compechem.systems import MDTrajectory
from compechem.tools import process_output
from compechem.tools import save_dftb_trajectory
from compechem.tools import compress_dftb_trajectory
import logging

logger = logging.getLogger(__name__)


class DFTBInput:
    """Interface for running DFTB+ calculations"""

    def __init__(
        self,
        hamiltonian: str = "DFTB",
        parameters: str = "3ob/3ob-3-1/",
        solver: str = None,
        dispersion: bool = False,
        parallel: str = "mpi",
        verbose: bool = False,
    ) -> None:
        """
        Parameters
        ----------
        hamiltonian : str, optional
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

        self.hamiltonian = hamiltonian
        self.parameters = parameters
        self.solver = solver
        self.dispersion = dispersion
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

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.hamiltonian} SPE")
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.hamiltonian.split()[0]}_spe",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            mol.write_gen(f"{mol.name}.gen")

            with open(f"{mol.name}.gen") as file:
                lines = file.readlines()
                atom_types = lines[1].split()

            with open("dftb_in.hsd", "w") as inp:

                inp.write(
                    "Geometry = GenFormat {\n"
                    f'  <<< "{mol.name}.gen"\n'
                    "}\n"
                    "\n"
                    "Driver = GeometryOptimization{\n"
                    "  MaxSteps = 0\n"
                    "}\n"
                    "\n"
                    f"Hamiltonian = {self.hamiltonian} {{\n"
                )

                if self.hamiltonian == "DFTB":
                    if self.solver:
                        inp.write(f"  Solver = {self.solver} {{}}\n")
                    inp.write(
                        "  Scc = Yes\n"
                        "  SlaterKosterFiles = Type2FileNames {\n"
                        f'    Prefix = "{self.parameters}"\n'
                        '    Separator = "-"\n'
                        '    Suffix = ".skf"\n'
                        "  }\n"
                        "  MaxAngularMomentum {\n"
                    )
                    for atom in atom_types:
                        inp.write(f'    {atom} = "{self.atom_dict[atom]}"\n')
                    inp.write("  }\n")
                    if mol.periodic:
                        inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    if "3ob" in self.parameters:
                        inp.write("  ThirdOrderFull = Yes\n" "  HubbardDerivs {\n")
                    for atom in atom_types:
                        inp.write(f"    {atom} = {self.hubbard_derivs[atom]}\n")
                    inp.write(
                        "  }\n"
                        "  HCorrection = Damping {\n"
                        "    Exponent = 4.00\n"
                        "  }\n"
                    )
                    if self.dispersion:
                        inp.write(
                            "  Dispersion = SimpleDftD3 {\n"
                            "    a1 = 0.746\n"
                            "    a2 = 4.191\n"
                            "    s6 = 1.0\n"
                            "    s8 = 3.209\n"
                            "  }\n"
                        )
                    inp.write("}\n")

                elif self.hamiltonian == "xTB":
                    if self.solver:
                        inp.write(f"  Solver = {self.solver} {{}}\n")
                    self.parameters = "gfn2"
                    inp.write('  Method = "GFN2-xTB"\n')
                    if mol.periodic:
                        inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    inp.write("}\n")

                inp.write("\n" "ParserOptions {\n" "  ParserVersion = 11\n" "}")

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > {self.output_path} 2>> output.err")

            elif self.parallel == "nompi":
                os.system(f"dftb+ > {self.output_path} 2>> output.err")

            with open("output.out", "r") as out:
                for line in out:
                    if "Total Energy" in line:
                        electronic_energy = float(line.split()[2])

            vibronic_energy = None

            if self.parameters in mol.energies:
                vibronic_energy = mol.energies[self.parameters].vibronic

            if inplace is False:

                mol.write_xyz(f"{mol.name}.xyz")

                newmol = System(f"{mol.name}.xyz", charge, spin, mol.periodic, mol.box_side)

                newmol.energies = copy.copy(mol.energies)

                newmol.energies[self.parameters] = Energies(
                    method=self.parameters,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

            else:
                mol.energies[self.parameters] = Energies(
                    method=self.parameters,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

            process_output(mol, self.hamiltonian, "spe", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def md_nvt(
        self,
        mol: System,
        steps: int,
        timestep: float = 1.0,
        temperature: int = 298,
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
        temperature : int, optional
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

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.hamiltonian} NVT MD")
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.hamiltonian.split()[0]}_md_nvt",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            mol.write_gen(f"{mol.name}.gen", box_side)

            with open(f"{mol.name}.gen") as file:
                lines = file.readlines()
                atom_types = lines[1].split()

            with open("dftb_in.hsd", "w") as inp:

                inp.write(
                    "Geometry = GenFormat {\n"
                    f'  <<< "{mol.name}.gen"\n'
                    "}\n"
                    "\n"
                    "Driver = VelocityVerlet{\n"
                    f"  TimeStep [fs] = {timestep}\n"
                    "  Thermostat = NoseHoover {\n"
                    f"    Temperature [Kelvin] = {temperature}\n"
                    "    CouplingStrength [cm^-1] = 3200\n"
                    "  }\n"
                    f"  Steps = {steps}\n"
                    "  MovedAtoms = 1:-1\n"
                    f"  MDRestartFrequency = {mdrestartfreq}\n"
                    "    Velocities [AA/ps] {\n"
                )
                for velocity in mol.velocities:
                    inp.write(f"      {velocity[1:]}")
                inp.write("  }\n" "}\n" "\n" f"Hamiltonian = {self.hamiltonian} {{\n")

                if self.hamiltonian == "DFTB":
                    if self.solver:
                        inp.write(f"  Solver = {self.solver} {{}}\n")
                    inp.write(
                        "  Scc = Yes\n"
                        "  SlaterKosterFiles = Type2FileNames {\n"
                        f'    Prefix = "{self.parameters}"\n'
                        '    Separator = "-"\n'
                        '    Suffix = ".skf"\n'
                        "  }\n"
                        "  MaxAngularMomentum {\n"
                    )
                    for atom in atom_types:
                        inp.write(f'    {atom} = "{self.atom_dict[atom]}"\n')
                    inp.write("  }\n")
                    if mol.periodic:
                        inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    if "3ob" in self.parameters:
                        inp.write("  ThirdOrderFull = Yes\n" "  HubbardDerivs {\n")
                    for atom in atom_types:
                        inp.write(f"    {atom} = {self.hubbard_derivs[atom]}\n")
                    inp.write(
                        "  }\n"
                        "  HCorrection = Damping {\n"
                        "    Exponent = 4.00\n"
                        "  }\n"
                    )
                    if self.dispersion:
                        inp.write(
                            "  Dispersion = SimpleDftD3 {\n"
                            "    a1 = 0.746\n"
                            "    a2 = 4.191\n"
                            "    s6 = 1.0\n"
                            "    s8 = 3.209\n"
                            "  }\n"
                        )
                    inp.write("}\n")

                elif self.hamiltonian == "xTB":
                    if self.solver:
                        inp.write(f"  Solver = {self.solver} {{}}\n")
                    self.parameters = "gfn2"
                    inp.write('  Method = "GFN2-xTB"\n')
                    if mol.periodic:
                        inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    inp.write("}\n")

                inp.write("\n" "ParserOptions {\n" "  ParserVersion = 11\n" "}")

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > {self.output_path} 2>> output.err")
            elif self.parallel == "nompi":
                os.environ["OMP_NUM_THREADS"] = f"{ncores}"
                os.system(f"dftb+ > {self.output_path} 2>> output.err")

            import random, string

            if inplace is False:
                newmol = System("geo_end.xyz", charge, spin, mol.periodic, mol.box_side)
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

            process_output(mol, self.hamiltonian, "md_nvt", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

        trajectory = MDTrajectory(f"{mol.name}_{suffix}", self.parameters)

        return trajectory

    def simulated_annealing(
        self,
        mol: System,
        start_temp: float = 100.0,
        target_temp: float = 5500.0,
        ramp_steps: int = 500,
        hold_steps: int = 200,
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
            Starting temperature (default, 100K)
        target_temp: float, optional
            Maximum temperature reached during the simulation (default, 5500K)
        ramp_steps: int, optional
            Number of MD steps for the heating/cooling ramps (default, 500 steps)
        hold_steps: int, optional
            Number of MD steps held at target_temp (default, 200 steps)
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
            f"{mol.name}, charge {charge} spin {spin} - {self.hamiltonian} Simulated Annealing"
        )
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.hamiltonian.split()[0]}_anneal",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            mol.write_gen(f"{mol.name}.gen", box_side)

            with open(f"{mol.name}.gen") as file:
                lines = file.readlines()
                atom_types = lines[1].split()

            with open("dftb_in.hsd", "w") as inp:

                inp.write(
                    "Geometry = GenFormat {\n"
                    f'  <<< "{mol.name}.gen"\n'
                    "}\n"
                    "\n"
                    "Driver = VelocityVerlet{\n"
                    f"  TimeStep [fs] = {timestep}\n"
                    "  Thermostat = NoseHoover {\n"
                    "    Temperature [Kelvin] = TemperatureProfile{\n"
                    f"      constant 1 {start_temp}\n"
                    f"      linear {ramp_steps-1} {target_temp}\n"
                    f"      constant {hold_steps} {target_temp}\n"
                    f"      linear {ramp_steps} {start_temp}\n"
                    "    }\n"
                    "    CouplingStrength [cm^-1] = 3200\n"
                    "  }\n"
                    "  MovedAtoms = 1:-1\n"
                    f"  MDRestartFrequency = {mdrestartfreq}\n"
                    "    Velocities [AA/ps] {\n"
                )
                for velocity in mol.velocities:
                    inp.write(f"      {velocity[1:]}")
                inp.write("  }\n" "}\n" "\n" f"Hamiltonian = {self.hamiltonian} {{\n")

                if self.hamiltonian == "DFTB":
                    if self.solver:
                        inp.write(f"  Solver = {self.solver} {{}}\n")
                    inp.write(
                        "  Scc = Yes\n"
                        "  SlaterKosterFiles = Type2FileNames {\n"
                        f'    Prefix = "{self.parameters}"\n'
                        '    Separator = "-"\n'
                        '    Suffix = ".skf"\n'
                        "  }\n"
                        "  MaxAngularMomentum {\n"
                    )
                    for atom in atom_types:
                        inp.write(f'    {atom} = "{self.atom_dict[atom]}"\n')
                    inp.write("  }\n")
                    if mol.periodic:
                        inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    if "3ob" in self.parameters:
                        inp.write("  ThirdOrderFull = Yes\n" "  HubbardDerivs {\n")
                    for atom in atom_types:
                        inp.write(f"    {atom} = {self.hubbard_derivs[atom]}\n")
                    inp.write(
                        "  }\n"
                        "  HCorrection = Damping {\n"
                        "    Exponent = 4.00\n"
                        "  }\n"
                    )
                    if self.dispersion:
                        inp.write(
                            "  Dispersion = SimpleDftD3 {\n"
                            "    a1 = 0.746\n"
                            "    a2 = 4.191\n"
                            "    s6 = 1.0\n"
                            "    s8 = 3.209\n"
                            "  }\n"
                        )
                    inp.write("}\n")

                elif self.hamiltonian == "xTB":
                    if self.solver:
                        inp.write(f"  Solver = {self.solver} {{}}\n")
                    self.parameters = "gfn2"
                    inp.write('  Method = "GFN2-xTB"\n')
                    if mol.periodic:
                        inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    inp.write("}\n")

                inp.write("\n" "ParserOptions {\n" "  ParserVersion = 11\n" "}")

            if self.parallel == "mpi":
                os.environ["OMP_NUM_THREADS"] = "1"
                os.system(f"mpirun -np {ncores} dftb+ > {self.output_path} 2>> output.err")
            elif self.parallel == "nompi":
                os.environ["OMP_NUM_THREADS"] = f"{ncores}"
                os.system(f"dftb+ > {self.output_path} 2>> output.err")

            import random, string

            if inplace is False:
                newmol = System("geo_end.xyz", charge, spin, mol.periodic, mol.box_side)
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

            process_output(mol, self.hamiltonian, "anneal", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

        trajectory = MDTrajectory(f"{mol.name}_{suffix}", self.parameters)

        return trajectory
