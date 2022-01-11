import os, shutil
from tempfile import mkdtemp


class OrcaInput:
    def __init__(
        self,
        method: str,
        basis_set: str = "def2-TZVP",
        aux_basis: str = "def2/J",
        nproc: int = 1,
        maxcore: int = 350,
        solvation: bool = False,
        solvent: str = "water",
        optionals: str = "",
    ) -> None:

        self.method = method
        self.basis_set = basis_set
        self.aux_basis = aux_basis
        self.nproc = nproc
        self.maxcore = maxcore
        self.solvation = solvation
        self.solvent = solvent
        self.optionals = optionals

    def spe(self, mol, remove_tdir=True):

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name} - {self.method} SPE")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_spe",
            dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz("geom.xyz")

        with open("input.inp", "w") as inp:
            inp.write(
                f"%pal nproc {self.nproc} end\n"
                f"%maxcore {self.maxcore}\n"
                f"! {self.method} {self.basis_set} {self.optionals}\n"
                f"! RIJCOSX {self.aux_basis}\n"
            )
            if self.solvation is True:
                inp.write(
                    "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                )
            inp.write(f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n")

        os.system("$ORCADIR/orca input.inp > output.out")

        with open("output.out", "r") as out:
            for line in out:
                if "FINAL SINGLE POINT ENERGY" in line:
                    electronic_energy = line.split()[-1]

        mol.energies[f"{self.method}"] = mol.Energies(
            method=f"{self.method}", electronic=electronic_energy,
        )

        if remove_tdir is True:
            shutil.rmtree(tdir)
        os.chdir(parent_dir)

        return

    def opt(self, mol, remove_tdir=True):

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name} - {self.method} OPT")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_opt",
            dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz("geom.xyz")

        with open("input.inp", "w") as inp:
            inp.write(
                f"%pal nproc {self.nproc} end\n"
                f"%maxcore {self.maxcore}\n"
                f"! {self.method} {self.basis_set} {self.optionals}\n"
                f"! RIJCOSX {self.aux_basis}\n"
            )
            if self.solvation is True:
                inp.write(
                    "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                )
            inp.write(f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n")

        os.system("$ORCADIR/orca input.inp > output.out")

        if remove_tdir is True:
            shutil.rmtree(tdir)
        os.chdir(parent_dir)

        return

    def freq(self, mol, remove_tdir=True):

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name} - {self.method} FREQ")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_freq",
            dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz("geom.xyz")

        with open("input.inp", "w") as inp:
            inp.write(
                f"%pal nproc {self.nproc} end\n"
                f"%maxcore {self.maxcore}\n"
                f"! {self.method} {self.basis_set} {self.optionals}\n"
                f"! RIJCOSX {self.aux_basis}\n"
            )
            if self.solvation is True:
                inp.write(
                    "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                )
            inp.write(f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n")

        os.system("$ORCADIR/orca input.inp > output.out")

        with open("output.out", "r") as out:
            for line in out:
                if "FINAL SINGLE POINT ENERGY" in line:
                    electronic_energy = line.split()[-1]
                if "G-E(el)" in line:
                    vibronic_energy = line.split()[-4]

        mol.energies[f"{self.method}"] = mol.Energies(
            method=f"{self.method}",
            electronic=electronic_energy,
            vibronic=vibronic_energy,
        )

        if remove_tdir is True:
            shutil.rmtree(tdir)
        os.chdir(parent_dir)

        return

    def nfreq(self, mol, remove_tdir=True):

        parent_dir = os.getcwd()
        print(f"INFO: {mol.name} - {self.method} NFREQ")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.method.split()[0]}_nfreq",
            dir=os.getcwd(),
        )

        os.chdir(tdir)
        mol.write_xyz("geom.xyz")

        with open("input.inp", "w") as inp:
            inp.write(
                f"%pal nproc {self.nproc} end\n"
                f"%maxcore {self.maxcore}\n"
                f"! {self.method} {self.basis_set} {self.optionals}\n"
                f"! RIJCOSX {self.aux_basis}\n"
            )
            if self.solvation is True:
                inp.write(
                    "%CPCM\n" "  SMD True\n" f'  SMDsolvent "{self.solvent}"\n' "end\n"
                )
            inp.write(f"* xyzfile {mol.charge} {mol.spin} geom.xyz\n")

        os.system("$ORCADIR/orca input.inp > output.out")

        with open("output.out", "r") as out:
            for line in out:
                if "FINAL SINGLE POINT ENERGY" in line:
                    electronic_energy = line.split()[-1]
                if "G-E(el)" in line:
                    vibronic_energy = line.split()[-4]

        mol.energies[f"{self.method}"] = mol.Energies(
            method=f"{self.method}",
            electronic=electronic_energy,
            vibronic=vibronic_energy,
        )

        if remove_tdir is True:
            shutil.rmtree(tdir)
        os.chdir(parent_dir)

        return

