import os, copy, shutil, sh
from tempfile import mkdtemp
from compechem.config import get_ncores
from compechem.systems import Ensemble, System, Energies
from compechem.tools import process_output
import logging

logger = logging.getLogger(__name__)


class OrcaInput:
    """Interface for running Orca calculations."""

    def __init__(
        self,
        method: str = "PBE",
        basis_set: str = "def2-TZVP",
        aux_basis: str = "def2/J",
        solvation: bool = False,
        solvent: str = "water",
        optionals: str = "",
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
        solvation : bool, optional
            CPCM(SMD) implicit solvation model, by default False
        solvent : str, optional
            SMD solvent, by default "water"
        optionals : str, optional
            optional keywords, by default ""
        """

        self.method = method
        self.basis_set = basis_set
        self.aux_basis = aux_basis
        self.solvation = solvation
        self.solvent = solvent
        self.optionals = optionals

    def spe(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
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

            mol.write_xyz(f"{mol.name}.xyz")

            with open("input.inp", "w") as inp:
                inp.write(
                    f"%pal nprocs {ncores} end\n"
                    f"%maxcore {maxcore}\n"
                    f"! {self.method} {self.basis_set} {self.optionals}\n"
                )
                if self.aux_basis:
                    inp.write(f"! RIJCOSX {self.aux_basis}\n")
                if self.solvation is True:
                    inp.write(
                        "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                    )
                inp.write(f"* xyzfile {charge} {spin} {mol.name}.xyz\n")

            os.system("$ORCADIR/orca input.inp > output.out")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])

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
        maxcore: int = 350,
        charge: int = None,
        spin: int = None,
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
            mol.write_xyz(f"{mol.name}.xyz")

            with open("input.inp", "w") as inp:
                inp.write(
                    f"%pal nprocs {ncores} end\n"
                    f"%maxcore {maxcore}\n"
                    f"! {self.method} {self.basis_set} {self.optionals}\n"
                )
                if self.aux_basis:
                    inp.write(f"! RIJCOSX {self.aux_basis}\n")
                if self.solvation is True:
                    inp.write(
                        "! Opt NumFreq\n"
                        "%CPCM\n"
                        "  SMD True\n"
                        f'  SMDsolvent "{self.solvent}"\n'
                        "end\n"
                    )
                else:
                    inp.write("! Opt Freq\n")
                inp.write(f"* xyzfile {charge} {spin} {mol.name}.xyz\n")

            os.system("$ORCADIR/orca input.inp > output.out")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])
                    if "G-E(el)" in line:
                        vibronic_energy = float(line.split()[-4])

            if inplace is False:

                newmol = System(f"{mol.name}.xyz", charge, spin)

                newmol.energies = copy.copy(mol.energies)

                newmol.energies[self.method] = Energies(
                    method=self.method,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

                newmol.update_geometry(f"{mol.name}.xyz")

            else:
                mol.energies[self.method] = Energies(
                    method=self.method,
                    electronic=electronic_energy,
                    vibronic=vibronic_energy,
                )

                mol.update_geometry(f"{mol.name}.xyz")

            process_output(mol, self.method, "opt", charge, spin)
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

            mol.write_xyz(f"{mol.name}.xyz")

            with open("input.inp", "w") as inp:
                inp.write(
                    f"%pal nprocs {ncores} end\n"
                    f"%maxcore {maxcore}\n"
                    f"! {self.method} {self.basis_set} {self.optionals}\n"
                )
                if self.aux_basis:
                    inp.write(f"! RIJCOSX {self.aux_basis}\n")
                if self.solvation is True:
                    inp.write(
                        "! NumFreq\n"
                        "%CPCM\n"
                        "  SMD True\n"
                        f'  SMDsolvent "{self.solvent}"\n'
                        "end\n"
                    )
                else:
                    inp.write("! Freq\n")
                inp.write(f"* xyzfile {charge} {spin} {mol.name}.xyz\n")

            os.system("$ORCADIR/orca input.inp > output.out")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])
                    if "G-E(el)" in line:
                        vibronic_energy = float(line.split()[-4])

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

            mol.write_xyz(f"{mol.name}.xyz")

            with open("input.inp", "w") as inp:
                inp.write(
                    f"%pal nprocs {ncores} end\n"
                    f"%maxcore {maxcore}\n"
                    f"! {self.method} {self.basis_set} {self.optionals}\n"
                    "! NumFreq\n"
                )
                if self.aux_basis:
                    inp.write(f"! RIJCOSX {self.aux_basis}\n")
                if self.solvation is True:
                    inp.write(
                        "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                    )

                inp.write(f"* xyzfile {charge} {spin} {mol.name}.xyz\n")

            os.system("$ORCADIR/orca input.inp > output.out")

            with open("output.out", "r") as out:
                for line in out:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        electronic_energy = float(line.split()[-1])
                    if "G-E(el)" in line:
                        vibronic_energy = float(line.split()[-4])

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
            mol.write_xyz(f"{mol.name}.xyz")

            with open("input.inp", "w") as inp:
                inp.write(
                    f"%pal nprocs {ncores} end\n"
                    f"%maxcore {maxcore}\n"
                    f"! {self.method} {self.basis_set} {self.optionals}\n"
                )
                if self.aux_basis:
                    inp.write(f"! RIJCOSX {self.aux_basis}\n")
                inp.write("! Opt\n")
                if self.solvation is True:
                    inp.write(
                        "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                    )
                inp.write("%geom\n" "  scan\n" f"    {scan}\n" "  end\n")
                if constraints:
                    inp.write("  constraints\n" f"    {{ {constraints} C }}\n" "  end\n")
                if invertconstraints:
                    inp.write("  invertConstraints true\n")
                inp.write("end\n")
                inp.write(f"* xyzfile {charge} {spin} {mol.name}.xyz\n")

            os.system("$ORCADIR/orca input.inp > output.out")

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
            solvation=True,
            solvent="water",
            optionals="",
        )


class r2SCAN(OrcaInput):
    def __init__(self):
        super().__init__(
            method="r2SCAN-3c",
            basis_set="",
            aux_basis=None,
            solvation=True,
            solvent="water",
            optionals="",
        )


class CCSD(OrcaInput):
    def __init__(self):
        super().__init__(
            method="DLPNO-CCSD",
            basis_set="Extrapolate(2/3,ANO)",
            aux_basis="AutoAux",
            solvation=True,
            solvent="water",
            optionals="",
        )
