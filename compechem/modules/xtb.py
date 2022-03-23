import os
from tempfile import mkdtemp
from compechem.modules import tools
from compechem.molecule import Molecule


class XtbInput:
    """Interface for running xTB calculations
    """

    def __init__(
        self,
        method: str = "gfn2",
        nproc: int = 1,
        solvation: bool = True,
        solvent: str = "water",
        optionals: str = "",
    ) -> None:
        """
        Parameters
        ----------
        method : str, optional
            level of theory, by default "gfn2"
        nproc : int, optional
            number of cores, by default 1
        solvation : bool, optional
            ALPB implicit solvation model, by default True
        solvent : str, optional
            ALPB solvent, by default "water"
        optionals : str, optional
            optional keywords/flags, by default ""
        """

        self.method = method

        self.nproc = nproc
        self.solvation = solvation
        self.solvent = solvent
        self.optionals = optionals

    def spe(self, mol, charge=None, spin=None, inplace=False, remove_tdir=True):
        """Single point energy calculation.

        Parameters
        ----------
        mol : Molecule object
            Input molecule to use in the calculation.
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

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name}, charge {charge} spin {spin} - {self.method} SPE")

        tdir = mkdtemp(
            prefix=mol.name + "_", suffix=f"_{self.method.split()[0]}_spe", dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz(f"{mol.name}.xyz")

        if self.solvation is True:
            os.system(
                f"xtb {mol.name}.xyz --{self.method} --alpb {self.solvent} --chrg {charge} --uhf {spin-1} -P {self.nproc} {self.optionals} > output.out 2>> output.err"
            )

        else:
            os.system(
                f"xtb {mol.name}.xyz --{self.method} --chrg {charge} --uhf {spin-1} -P {self.nproc} {self.optionals} > output.out 2>> output.err"
            )

        with open("output.out", "r") as out:
            for line in out:
                if "TOTAL ENERGY" in line:
                    electronic_energy = float(line.split()[-3])

        vibronic_energy = None

        if self.method in mol.energies:
            vibronic_energy = mol.energies[f"{self.method}"].vibronic

        if inplace is False:

            newmol = Molecule(f"{mol.name}.xyz", charge, spin)

            newmol.energies = mol.energies

            newmol.energies[f"{self.method}"] = newmol.Energies(
                method=f"{self.method}", electronic=electronic_energy, vibronic=vibronic_energy
            )

        else:
            mol.energies[f"{self.method}"] = mol.Energies(
                method=f"{self.method}", electronic=electronic_energy, vibronic=vibronic_energy
            )

        tools.process_output(mol, self.method, charge, spin, "spe", tdir, remove_tdir, parent_dir)

        if inplace is False:
            return newmol

    def opt(self, mol, charge=None, spin=None, inplace=False, remove_tdir=True):
        """Geometry optimization + frequency analysis.

        Parameters
        ----------
        mol : Molecule object
            Input molecule to use in the calculation
        charge : int, optional
            Total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            Total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            Temporary work directory will be removed, by default True

        Returns
        -------
        newmol : Molecule object
            Output molecule containing the new geometry and energies.
        
        If a dissociation or a cyclization is observed, ignore the calculation and return the 
        original "mol" molecule.
        """

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name}, charge {charge} spin {spin} - {self.method} OPT")

        tdir = mkdtemp(
            prefix=mol.name + "_", suffix=f"_{self.method.split()[0]}_opt", dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz(f"{mol.name}.xyz")

        if self.solvation is True:
            os.system(
                f"xtb {mol.name}.xyz --{self.method} --alpb {self.solvent} --chrg {charge} --uhf {spin-1} --ohess -P {self.nproc} {self.optionals} > output.out 2>> output.err"
            )

        else:
            os.system(
                f"xtb {mol.name}.xyz --{self.method} --chrg {charge} --uhf {spin-1} --ohess -P {self.nproc} {self.optionals} > output.out 2>> output.err"
            )

        try: 
            tools.dissociation_check(mol)
        except:
            print(f"ERROR: Exception occcurred for {mol.name}. Ignoring optimization.")
            return mol

        tools.cyclization_check(mol, f"{mol.name}.xyz", "xtbopt.xyz")

        with open("output.out", "r") as out:
            for line in out:
                if "TOTAL ENERGY" in line:
                    electronic_energy = float(line.split()[-3])
                if "G(RRHO) contrib." in line:
                    vibronic_energy = float(line.split()[-3])

        if inplace is False:

            newmol = Molecule(f"{mol.name}.xyz", charge, spin)

            newmol.energies = mol.energies

            newmol.energies[f"{self.method}"] = newmol.Energies(
                method=f"{self.method}", electronic=electronic_energy, vibronic=vibronic_energy
            )

            newmol.update_geometry("xtbopt.xyz")

        else:
            mol.energies[f"{self.method}"] = mol.Energies(
                method=f"{self.method}", electronic=electronic_energy, vibronic=vibronic_energy
            )

            mol.update_geometry("xtbopt.xyz")

        tools.process_output(mol, self.method, charge, spin, "opt", tdir, remove_tdir, parent_dir)

        if inplace is False:
            return newmol

    def freq(self, mol, charge=None, spin=None, inplace=False, remove_tdir=True):
        """Frequency analysis.

        Parameters
        ----------
        mol : Molecule object
            input molecule to use in the calculation
        charge : int, optional
            Total charge of the molecule. Default is taken from the input molecule.
        spin : int, optional
            Total spin of the molecule. Default is taken from the input molecule.
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True

        Returns
        -------
        newmol : Molecule object
            Output molecule containing the new energies.
        """

        if charge is None:
            charge = mol.charge
        if spin is None:
            spin = mol.spin

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name}, charge {charge} spin {spin} - {self.method} FREQ")

        tdir = mkdtemp(
            prefix=mol.name + "_", suffix=f"_{self.method.split()[0]}_freq", dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz(f"{mol.name}.xyz")

        if self.solvation is True:
            os.system(
                f"xtb {mol.name}.xyz --{self.method} --alpb {self.solvent} --chrg {charge} --uhf {spin-1} --hess -P {self.nproc} {self.optionals} > output.out 2>> output.err"
            )

        else:
            os.system(
                f"xtb {mol.name}.xyz --{self.method} --chrg {charge} --uhf {spin-1} --hess -P {self.nproc} {self.optionals} > output.out 2>> output.err"
            )

        with open("output.out", "r") as out:
            for line in out:
                if "TOTAL ENERGY" in line:
                    electronic_energy = float(line.split()[-3])
                if "G(RRHO) contrib." in line:
                    vibronic_energy = float(line.split()[-3])

        if inplace is False:

            newmol = Molecule(f"{mol.name}.xyz", charge, spin)

            newmol.energies = mol.energies

            newmol.energies[f"{self.method}"] = newmol.Energies(
                method=f"{self.method}", electronic=electronic_energy, vibronic=vibronic_energy
            )

        else:
            mol.energies[f"{self.method}"] = mol.Energies(
                method=f"{self.method}", electronic=electronic_energy, vibronic=vibronic_energy
            )

        tools.process_output(mol, self.method, charge, spin, "freq", tdir, remove_tdir, parent_dir)

        if inplace is False:
            return newmol
