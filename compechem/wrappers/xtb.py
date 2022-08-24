import os, copy, sh, shutil
from tempfile import mkdtemp
from compechem.config import get_ncores
from compechem.systems import System, Energies
from compechem.tools import process_output
from compechem.tools import add_flag
from compechem.tools import dissociation_check
from compechem.tools import cyclization_check
import logging

logger = logging.getLogger(__name__)


class XtbInput:
    """Interface for running xTB calculations"""

    def __init__(
        self,
        method: str = "gfn2",
        solvation: bool = True,
        solvent: str = "water",
        optionals: str = "",
    ) -> None:
        """
        Parameters
        ----------
        method : str, optional
            level of theory, by default "gfn2"
        solvation : bool, optional
            ALPB implicit solvation model, by default True
        solvent : str, optional
            ALPB solvent, by default "water"
        optionals : str, optional
            optional keywords/flags, by default ""
        """

        self.method = method
        self.solvation = solvation
        self.solvent = solvent
        self.optionals = optionals

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
            Input system to use in the calculation.
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            total charge of the system. Default is taken from the input system.
        spin : int, optional
            total spin of the system. Default is taken from the input system.
        inplace : bool, optional
            updates info for the input system instead of outputting a new system object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output system containing the new energies.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} SPE")
        logger.debug(f"Running xTB calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_spe",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            mol.write_xyz(f"{mol.name}.xyz")

            if self.solvation is True:
                os.system(
                    f"xtb {mol.name}.xyz --{self.method} --alpb {self.solvent} --chrg {charge} --uhf {spin-1} -P {ncores} {self.optionals} > output.out 2>> output.err"
                )

            else:
                os.system(
                    f"xtb {mol.name}.xyz --{self.method} --chrg {charge} --uhf {spin-1} -P {ncores} {self.optionals} > output.out 2>> output.err"
                )

            with open("output.out", "r") as out:
                for line in out:
                    if "TOTAL ENERGY" in line:
                        electronic_energy = float(line.split()[-3])

            vibronic_energy = None

            if self.method in mol.energies:
                vibronic_energy = mol.energies[self.method].vibronic

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

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
        ncores: int = None,
        maxcore=None,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Geometry optimization + frequency analysis.

        Parameters
        ----------
        mol : System object
            Input system to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            Total charge of the system. Default is taken from the input system.
        spin : int, optional
            Total spin of the system. Default is taken from the input system.
        inplace : bool, optional
            updates info for the input system instead of outputting a new system object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output system containing the new geometry and energies.

        If a dissociation or a cyclization is observed, ignore the calculation and return the
        original "mol" system.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} OPT")
        logger.debug(f"Running xTB calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_opt",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            mol.write_xyz(f"{mol.name}.xyz")

            if self.solvation is True:
                os.system(
                    f"xtb {mol.name}.xyz --{self.method} --alpb {self.solvent} --chrg {charge} --uhf {spin-1} --ohess -P {ncores} {self.optionals} > output.out 2>> output.err"
                )

            else:
                os.system(
                    f"xtb {mol.name}.xyz --{self.method} --chrg {charge} --uhf {spin-1} --ohess -P {ncores} {self.optionals} > output.out 2>> output.err"
                )

            if dissociation_check() is True:
                logger.error(f"Dissociation spotted for {mol.name}.")
                add_flag(
                    mol,
                    f"Dissociation occurred during geometry optimization with {self.method}.",
                )
                return None
            elif cyclization_check(f"{mol.name}.xyz", "xtbopt.xyz") is True:
                logger.error(f"Cyclization change spotted for {mol.name}.")
                add_flag(
                    mol,
                    f"Cyclization change occurred during geometry optimization with {self.method}.",
                )
                return None
            else:
                with open("output.out", "r") as out:
                    for line in out:
                        if "TOTAL ENERGY" in line:
                            electronic_energy = float(line.split()[-3])
                        if "G(RRHO) contrib." in line:
                            vibronic_energy = float(line.split()[-3])

                if inplace is False:

                    newmol = System(f"{mol.name}.xyz", charge, spin)

                    newmol.energies = copy.copy(mol.energies)

                    newmol.energies[self.method] = Energies(
                        method=self.method,
                        electronic=electronic_energy,
                        vibronic=vibronic_energy,
                    )

                    newmol.update_geometry("xtbopt.xyz")

                else:
                    mol.energies[self.method] = Energies(
                        method=self.method,
                        electronic=electronic_energy,
                        vibronic=vibronic_energy,
                    )

                    mol.update_geometry("xtbopt.xyz")

                process_output(mol, self.method, "opt", charge, spin)
                if remove_tdir:
                    shutil.rmtree(tdir)

        if inplace is False:
            return newmol

    def freq(
        self,
        mol: System,
        ncores: int = None,
        maxcore=None,
        charge: int = None,
        spin: int = None,
        inplace: bool = False,
        remove_tdir: bool = True,
    ):
        """Frequency analysis.

        Parameters
        ----------
        mol : System object
            input system to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : dummy variable
            dummy variable used for compatibility with Orca calculations
        charge : int, optional
            Total charge of the system. Default is taken from the input system.
        spin : int, optional
            Total spin of the system. Default is taken from the input system.
        inplace : bool, optional
            updates info for the input system instead of outputting a new system object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        newmol : System object
            Output system containing the new energies.
        """

        if ncores is None:
            ncores = get_ncores()

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        logger.info(f"{mol.name}, charge {charge} spin {spin} - {self.method} FREQ")
        logger.debug(f"Running xTB calculation on {ncores} cores")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_freq",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            mol.write_xyz(f"{mol.name}.xyz")

            if self.solvation is True:
                os.system(
                    f"xtb {mol.name}.xyz --{self.method} --alpb {self.solvent} --chrg {charge} --uhf {spin-1} --hess -P {ncores} {self.optionals} > output.out 2>> output.err"
                )

            else:
                os.system(
                    f"xtb {mol.name}.xyz --{self.method} --chrg {charge} --uhf {spin-1} --hess -P {ncores} {self.optionals} > output.out 2>> output.err"
                )

            with open("output.out", "r") as out:
                for line in out:
                    if "TOTAL ENERGY" in line:
                        electronic_energy = float(line.split()[-3])
                    if "G(RRHO) contrib." in line:
                        vibronic_energy = float(line.split()[-3])

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

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

            process_output(mol, self.method, "freq", charge, spin)
            if remove_tdir:
                shutil.rmtree(tdir)

        if inplace is False:
            return newmol
