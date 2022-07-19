import os, copy, sys
from tempfile import mkdtemp
from compechem.config import get_ncores
from compechem.molecule import Molecule, Energies
from compechem import tools
import logging

logger = logging.getLogger(__name__)


class DFTBInput:
    """Interface for running DFTB+ calculations"""

    def __init__(
        self,
        geom_type: str = "C",
        box_side: float = None,
        hamiltonian: str = "DFTB",
        parameters: str = "3ob/3ob-3-1/",
        dispersion: bool = False,
        parallel: str = "mpi",
    ) -> None:
        """
        Parameters
        ----------
        geom_type : str
            type of geometry. C = cluster (single molecule), S = supercell (periodic system)
        box_side : float
            size of the periodic box size (in Angstrom).
        hamiltonian : str, optional
            level of theory, by default "DFTB". "xTB" also supported.
        parameters : str, optional
            parameters to be used for the DFTB Hamiltonian (by default 3ob)
        dispersion : bool, optional
            activates D3 dispersion corrections (off by default)
        parallel : str, optional
            selects either openmpi-parallel version (mpi) or shared memory version (nompi)
        """

        self.hamiltonian = hamiltonian
        self.parameters = parameters
        self.dispersion = dispersion
        self.geom_type = geom_type
        self.box_side = box_side
        self.parallel = parallel

        self.atom_dict = {
            "H": "s",
            "O": "p",
            "C": "p",
            "N": "p",
            "S": "d",
            "F": "p",
            "Cl": "d",
            "Br": "d",
            "I": "d",
        }

        self.hubbard_derivs = {
            "C": -0.1492,
            "O": -0.1575,
            "H": -0.1857,
            "N": -0.1535,
        }

    def spe(
        self,
        mol: Molecule,
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
        mol : Molecule object
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
        newmol : Molecule object
            Output molecule containing the new energies.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        parent_dir = os.getcwd()
        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.hamiltonian} SPE")
        logger.debug(f"Running DFTB+ calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.hamiltonian.split()[0]}_spe",
            dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_gen(f"{mol.name}.gen", geom_type=self.geom_type, box_side=self.box_side)

        with open(f"{mol.name}.gen") as file:
            lines = file.readlines()
            atom_types = lines[1].split()

        with open("dftb_in.hsd", "w") as inp:

            inp.write(
                f"""Geometry = GenFormat {{
  <<< "{mol.name}.gen"
}}

Driver = GeometryOptimization{{
  MaxSteps = 0
}}

Hamiltonian = {self.hamiltonian} {{"""
            )
            if self.hamiltonian == "DFTB":
                inp.write(
                    f"""
  Scc = Yes
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{self.parameters}"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
"""
                )
            for atom in atom_types:
                inp.write(f'    {atom} = "{self.atom_dict[atom]}"\n')
            inp.write("  }\n")
            if self.geom_type == "S":
                inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
            if self.parameters == "3ob/3ob-3-1/":
                inp.write(
                    """  ThirdOrderFull = Yes
  HubbardDerivs {
"""
                )
                for atom in atom_types:
                    inp.write(f"    {atom} = {self.hubbard_derivs[atom]}\n")
                inp.write(
                    """  }
  HCorrection = Damping {
    Exponent = 4.00
  }
"""
                )
            if self.dispersion:
                inp.write(
                    """  Dispersion = SimpleDftD3 {
    a1 = 0.746
    a2 = 4.191
    s6 = 1.0
    s8 = 3.209
  }"""
                )
            inp.write("}\n")
            if self.hamiltonian == "xTB":
                self.parameters = "gfn2"
                inp.write('  Method = "GFN2-xTB"\n')
                if self.geom_type == "S":
                    inp.write("  kPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n")
                    inp.write("}\n")
            inp.write(
                """
ParserOptions {
  ParserVersion = 11
}"""
            )

        exit()

        if self.parallel == "mpi":
            os.system(f"mpirun -np {ncores} dftb+ > output.out 2>> output.err")
        elif self.parallel == "nompi":
            os.system(f"dftb+ > output.out 2>> output.err")

        with open("output.out", "r") as out:
            for line in out:
                if "Total Energy" in line:
                    print(line.split())
                    electronic_energy = float(line.split()[2])
                    print(electronic_energy)

        vibronic_energy = None

        if self.parameters in mol.energies:
            vibronic_energy = mol.energies[f"{self.parameters}"].vibronic

        ### NEEDS TO BE FIXED vvv
        if inplace is False:

            mol.write_xyz(f"{mol.name}.xyz")

            newmol = Molecule(f"{mol.name}.xyz", charge, spin)

            newmol.energies = copy.copy(mol.energies)

            newmol.energies[f"{self.parameters}"] = Energies(
                method=f"{self.parameters}",
                electronic=electronic_energy,
                vibronic=vibronic_energy,
            )
        ### NEEDS TO BE FIXED ^^^

        else:
            mol.energies[f"{self.parameters}"] = Energies(
                method=f"{self.parameters}",
                electronic=electronic_energy,
                vibronic=vibronic_energy,
            )

        tools.process_output(
            mol, self.hamiltonian, charge, spin, "spe", tdir, remove_tdir, parent_dir
        )

        if inplace is False:
            return newmol
