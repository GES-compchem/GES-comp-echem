import os, copy, shutil, sh, logging
import numpy as np
from typing import Dict, Optional, Any
from tempfile import mkdtemp

import compechem.config as cfg
from compechem.config import get_ncores
from compechem.systems import Ensemble, System
from compechem.tools import process_output
from compechem.tools.externalutilities import split_multixyz
from compechem.core.base import Engine
from compechem.core.dependency_finder import locate_orca
from compechem.core.spectroscopy import VibrationalData
from compechem.tools.internaltools import clean_suffix


logger = logging.getLogger(__name__)


class OrcaJobInfo:
    """
    The OrcaJobInfo class is an helper class designed to provide a compact tool to define orca-based quantum chemical
    calculations. The attributes of the class univocally identify a calculation type and can be used to compile orca
    input files.

    Attributes
    ----------
    is_singlet: bool
        If True, the molecule on which the calculation is performed is in a singlet state (If cube is saved also the 
        spindensity will be saved).
    solvent: Optional[str]
        The SMD solvent to be used during the calculation.
    opt: bool
        If set to True, will trigger a geometry optimization calculation.
    opt_ts: bool
        If set to True, will trigger a transition-state optimization calculation.
    freq: bool
        If set to True, will trigger a frequency analysis (with analytical frequencies if not in solvent).
    nfreq: bool
        If set to True, will trigger a numerical frequency analysis.
    scan: Optional[str]
        If not None, will encode the scan parameters in string format.
    scan_ts: Optional[str]
        If not None, will encode the transition state scan parameters in string format.
    neb_ci: bool
        If set to True will trigger a NEB-CI calculation.
    constraints: Optional[str]
        If not None, will encode the constraints to be used during the scan or optimization operations.
    invert_constraints: bool
        If set to True, will invert the constraints to be used in the calculation.
    fullscan: bool
        If set to True during scan TS will continue the scan even if the transition state is located.
    cube_dim: Optional[int]
        If set to integer value will trigger the saving of cube files with a equally spaced grid of cube_dim points.
    calc_hess: bool
        If set to True will compute the analytical hessian before the starting of a calculation.
    hirshfeld: bool
        If set to True will trigger the Hirshfeld population analysis
    nearir: bool
        If set to True will trigger the calcuation of overtones and combination bands in the infrared spectrum
    raman: bool
        If set to True, will trigger the polarizability calcualtion outputting also the raman spectra.
    neb_product: Optional[str]
        If not None, will encode the final endopoint to be used in the neb calculation.
    neb_images: Optional[int]
        If set to an integer value, will define the number of images to be used in the neb calculation (endpoints excluded)
    neb_preopt: bool
        If set to True when a NEB calculation is performed, will preoptimize the reactant and product structures.
    """
    def __init__(self) -> None:
        self.__ncores: int = get_ncores()
        self.__maxcore: int = 750
        self.is_singlet: bool = True
        self.solvent: Optional[str] = None
        self.opt: bool = False
        self.opt_ts: bool = False
        self.freq: bool = False
        self.nfreq: bool = False
        self.scan: Optional[str] = None
        self.scan_ts: Optional[str] = None
        self.neb_ci: bool = False

        self.constraints: Optional[str] = None
        self.invert_constraints: bool = False
        self.fullscan: bool = False

        self.cube_dim: Optional[int] = None

        self.calc_hess: bool = False
        self.hirshfeld: bool = False
        self.nearir: bool = False
        self.raman: bool = False

        self.neb_product: Optional[str] = None
        self.neb_images: Optional[int] = None
        self.neb_preopt: bool = False

        self.__user_blocks: Dict[str, Dict[str, Any]] = {}

        self.__print_level: Optional[str] = None
        self.__optimization_level: Optional[str] = None
        self.__scf_convergence_level: Optional[str] = None
        self.__scf_convergence_strategy: Optional[str] = None


    @property
    def user_blocks(self) -> Dict[str, Dict[str, Any]]:
        return self.__user_blocks

    @user_blocks.setter
    def user_blocks(self, value: Dict[str, Dict[str, Any]]):

        blocks = {}
        
        # set the user block variable imposing lowercase block-name and keys
        for name, block in value.items():
            blocks[name.lower()] = {}
            for key, value in block.items():
                blocks[name.lower()][key.lower()] = value

        self.__user_blocks = blocks

    @property
    def geom_block(self) -> Dict[str, Any]:
        """
        Returns the %geom block generated by merging the user provided geometry options with the ones specified by the
        library flags.

        Returns
        -------
        Dict[str, Any]
            The content of the geom block encoded in dictionary format.
        """
        block = {} if "geom" not in self.user_blocks else self.user_blocks["geom"]

        if self.calc_hess is True:
            block["calc_hess"] = "true"

        if self.fullscan is True:
            block["fullscan"] = "true"

        if self.scan is not None:
            block["scan"] = f"{self.scan}  end"

        if self.scan_ts is not None:
            block["scan"] = f"{self.scan_ts}  end"

        if self.constraints is not None:
            block["constraints"] = f"{{ {self.constraints} C }}  end"

        if self.invert_constraints is True:
            block["invertconstraints"] = "true"
        
        return block
    

    @property
    def cpcm_block(self) -> Dict[str, Any]:
        """
        Returns the %cpcm block generated by merging the user provided solvent options with the ones specified by the
        library flags.

        Returns
        -------
        Dict[str, Any]
            The content of the cpcm block encoded in dictionary format.
        """
        block = {} if "cpcm" not in self.user_blocks else self.user_blocks["cpcm"]

        if self.solvent:
            block["smd"] = "true"
            block["smdsolvent"] = f'"{self.solvent}"'
        
        return block
    

    @property
    def plots_block(self) -> Dict[str, Any]:
        """
        Returns the %plots block generated by merging the user provided plot options with the ones specified by the
        library flags.

        Returns
        -------
        Dict[str, Any]
            The content of the geom block encoded in dictionary format.
        """
        block = {} if "plots" not in self.user_blocks else self.user_blocks["plots"]

        if self.cube_dim is not None:
            block["format"] = "gaussian_cube"
            block["dim1"] = f"{self.cube_dim}"
            block["dim2"] = f"{self.cube_dim}"
            block["dim3"] = f"{self.cube_dim}"
            block['eldens("eldens.cube");'] = ""

            if self.is_singlet is False:
                block['spindens("spindens.cube");'] = ""
        
        return block
    
    @property
    def output_block(self) -> Dict[str, Any]:
        """
        Returns the %output block generated by merging the user provided output options with the ones specified by the
        library flags.

        Returns
        -------
        Dict[str, Any]
            The content of the output block encoded in dictionary format.
        """
        block = {} if "output" not in self.user_blocks else self.user_blocks["output"]

        if self.hirshfeld is True:
            block["print[p_hirshfeld]"] = "1"

        return block
       

    @property
    def elprop_block(self) -> Dict[str, Any]:
        """
        Returns the %elprop block generated by merging the user provided electronic properties options with the ones 
        specified by the library flags.

        Returns
        -------
        Dict[str, Any]
            The content of the elprop block encoded in dictionary format.
        """

        block = {} if "elprop" not in self.user_blocks else self.user_blocks["elprop"]

        if self.raman is True:
            block["polar"] = "1"

        return block
    
    @property
    def neb_block(self) -> Dict[str, Any]:
        """
        Returns the %neb block generated by merging the user provided neb options with the ones specified by the library
        flags.

        Returns
        -------
        Dict[str, Any]
            The content of the neb block encoded in dictionary format.
        """

        block = {} if "neb" not in self.user_blocks else self.user_blocks["neb"]

        if self.neb_ci is True:

            block["product"] = f'"{self.neb_product}"'

            if self.neb_images is not None:
                block["nimages"] = self.neb_images
            
            if self.neb_preopt is True:
                block["preopt"] = "true"

        return block
        
    @property
    def parsed_blocks(self) -> Dict[str, Dict[str, Any]]:
        """
        Returns the dictionary encoding all the blocks to be used in the definition of an orca input file. The function 
        merges the user provided options with the ones specified by the library flags.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding all the blocks to be used in the orca input file.
        """
        blocks = {}

        if self.geom_block != {}:
            blocks["geom"] = self.geom_block
        
        if self.cpcm_block != {}:
            blocks["cpcm"] = self.cpcm_block
        
        if self.plots_block != {}:
            blocks["plots"] = self.plots_block
        
        if self.output_block != {}:
            blocks["output"] = self.output_block

        if self.elprop_block != {}:
            blocks["elprop"] = self.elprop_block
        
        if self.neb_block != {}:
            blocks["neb"] = self.neb_block
        
        for key, value in self.user_blocks.items():
            
            if key not in ["geom", "cpcm", "plots", "output", "elprop", "neb"]:
                blocks[key] = value
        
        return blocks
    
        
    @property
    def ncores(self) -> int:
        """
        The number of cores to be used in the calculation

        Returns
        -------
        int
            The total number of cores
        """
        return self.__ncores

    @ncores.setter
    def ncores(self, value: Optional[int]) -> None:
        self.__ncores = get_ncores() if value is None else value

    @property
    def maxcore(self) -> int:
        """
        The amount of RAM per core to be used in the calculation

        Returns
        -------
        int
            The maxcore value to be used by orca
        """
        return self.__maxcore

    @maxcore.setter
    def maxcore(self, value: Optional[int]) -> None:
        self.__maxcore = 750 if value is None else value

    @property
    def print_level(self) -> Optional[str]:
        """
        The print level to be used during the calculation. The possible values are "MINIPRINT", "SMALLPRINT",
        "NORMALPRINT", "LARGEPRINT".

        Returns
        -------
        str
            The string encoding the print level to be used in the output.
        """
        return self.__print_level

    @print_level.setter
    def print_level(self, value: Optional[str]) -> None:
        if value is not None:
            ACCEPTED_VALUES = [
                "MINIPRINT",
                "SMALLPRINT",
                "NORMALPRINT",
                "LARGEPRINT",
            ]
            if value.upper() not in ACCEPTED_VALUES:
                raise ValueError(
                    f"`{value}` is not a valid optimization level. Must be one of {', '.join(ACCEPTED_VALUES)}"
                )

            self.__print_level = value.upper()

        else:
            self.__print_level = None

    @property
    def optimization_level(self) -> Optional[str]:
        """
        The optimization level to be used during the calculation. The possible values are "VERYTIGHTOPT", "TIGHTOPT",
        "NORMALOPT", "LOOSEOPT".

        Returns
        -------
        str
            The string encoding the optimization level to be used in the output.
        """
        return self.__optimization_level

    @optimization_level.setter
    def optimization_level(self, value: str) -> None:
        if value is not None:
            ACCEPTED_VALUES = [
                "VERYTIGHTOPT",
                "TIGHTOPT",
                "NORMALOPT",
                "LOOSEOPT",
            ]
            if value.upper() not in ACCEPTED_VALUES:
                raise ValueError(
                    f"`{value}` is not a valid optimization level. Must be one of {', '.join(ACCEPTED_VALUES)}"
                )

            self.__optimization_level = value.upper()

        else:
            self.__print_level = None

    @property
    def scf_convergence_level(self) -> Optional[str]:
        """
        The scf convergence level to be used during the calculation. The possible values are "NORMALSCF", "LOOSESCF",
        "SLOPPYSCF", "STRONGSCF", "TIGHTSCF", "VERYTIGHTSCF", "EXTREMESCF".

        Returns
        -------
        str
            The string encoding the scf convergence level to be used during the calculation.
        """
        return self.__scf_convergence_level

    @scf_convergence_level.setter
    def scf_convergence_level(self, value: str) -> None:
        if value is not None:
            ACCEPTED_VALUES = [
                "NORMALSCF",
                "LOOSESCF",
                "SLOPPYSCF",
                "STRONGSCF",
                "TIGHTSCF",
                "VERYTIGHTSCF",
                "EXTREMESCF",
            ]
            if value.upper() not in ACCEPTED_VALUES:
                raise ValueError(
                    f"`{value}` is not a valid SCF convergence level. Must be one of {', '.join(ACCEPTED_VALUES)}"
                )

            self.__scf_convergence_level = value.upper()

        else:
            self.__print_level = None

    @property
    def scf_convergence_strategy(self) -> Optional[str]:
        """
        The scf convergence strategy to be used during the calculation. The possible values are "EASYCONV", "NORMALCONV",
        "SLOWCONV", "VERYSLOWCONV", "FORCECONV", "IGNORECONV".

        Returns
        -------
        str
            The string encoding the scf convergence strategy to be used during the calculation.
        """
        return self.__scf_convergence_strategy

    @scf_convergence_strategy.setter
    def scf_convergence_strategy(self, value: str) -> None:
        if value is not None:
            ACCEPTED_VALUES = [
                "EASYCONV",
                "NORMALCONV",
                "SLOWCONV",
                "VERYSLOWCONV",
                "FORCECONV",
                "IGNORECONV",
            ]
            if value.upper() not in ACCEPTED_VALUES:
                raise ValueError(
                    f"`{value}` is not a valid convergence strategy. Must be one of {', '.join(ACCEPTED_VALUES)}"
                )

            self.__scf_convergence_strategy = value.upper()

        else:
            self.__print_level = None


class OrcaInput(Engine):
    """Interface for running Orca calculations.

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
    blocks : Dict[str, Dict[str, Any]]
        The dictionary of dictionaries encoding a series of custom blocks defined by the user
    ORCADIR: str, optional
        the path or environment variable containing the path to the ORCA folder. If set
        to None (default) the orca executable will be loaded automatically.
    """

    def __init__(
        self,
        method: str = "PBE",
        basis_set: str = "def2-TZVP",
        aux_basis: str = "def2/J",
        solvent: str = None,
        optionals: str = "",
        blocks: Dict[str, Dict[str, Any]] = {},
        ORCADIR: str = None,
    ) -> None:
        super().__init__(method)

        self.basis_set = basis_set if basis_set else ""
        self.aux_basis = aux_basis if aux_basis else ""
        self.solvent = solvent
        self.optionals = optionals
        self.blocks = blocks
        self.__ORCADIR = ORCADIR if ORCADIR else locate_orca(get_folder=True)

        self.level_of_theory += f""" | basis: {basis_set} | solvent: {solvent}"""

        self.__output_suffix = f"orca_{method}"
        self.__output_suffix += f"_{basis_set}" if basis_set else ""
        self.__output_suffix += f"_{solvent}" if solvent else "_vacuum"
        self.__output_suffix = clean_suffix(self.__output_suffix)

    def write_input(
        self,
        mol: System,
        job_info: OrcaJobInfo,
    ) -> None:
        mol.geometry.write_xyz(f"{mol.name}.xyz")

        logger.debug(f"Running ORCA calculation on {job_info.ncores} cores and {job_info.maxcore} MB of RAM")

        input = (
            "%pal\n"
            f"  nprocs {job_info.ncores}\n"
            "end\n\n"
            f"%maxcore {job_info.maxcore}\n\n"
            f"! {self.method} {self.basis_set} {self.optionals}\n"
        )

        if job_info.scf_convergence_strategy is not None or job_info.scf_convergence_level is not None:
            input += "! "

            if job_info.scf_convergence_level is not None:
                input += job_info.scf_convergence_level
                input += " "

            if job_info.scf_convergence_strategy is not None:
                input += job_info.scf_convergence_strategy
                input += " "
            input += "\n"

        if self.aux_basis:
            input += f"! RIJCOSX {self.aux_basis}\n\n"

        if job_info.print_level is not None:
            input += f"! {job_info.print_level}\n\n"

        if (job_info.opt is True or job_info.opt_ts is True) and job_info.freq is True and self.solvent is not None:
            logger.warning("Optimization with frequency in solvent was requested. Switching to numerical frequencies.")
            job_info.freq = False
            job_info.nfreq = True

        if job_info.opt is True:
            input += "! Opt\n" if job_info.optimization_level is None else f"! {job_info.optimization_level}\n"

        if job_info.scan is not None:
            input += "! Opt\n"

        if job_info.opt_ts is True:
            input += "! OptTS\n"

        if job_info.scan_ts is not None:
            input += "! ScanTS\n"

        if job_info.freq is True:
            if self.solvent:
                logger.warning("Analytical frequencies are not supported for the SMD solvent model.")

            input += "! Freq\n"

        if job_info.nfreq is True:
            input += "! NumFreq\n"

        if job_info.nearir is True:
            input += "! NearIR\n"
        
        if job_info.neb_ci is True:
            input += "! NEB-CI\n"

        if job_info.parsed_blocks != {}:

            for name, block in job_info.parsed_blocks.items():

                input += f"%{name}\n"

                for key, entry in block.items():
                    if type(entry) is dict:
                        input += f"  {key}\n"
                        for subkey, subvalue in entry.items():
                            input += f"    {subkey} {subvalue}\n"
                        input += "  end\n"
                    else:
                        input += f"  {key} {entry}\n"

                input += "end\n\n"
        
        input += f"* xyzfile {mol.charge} {mol.spin} {mol.name}.xyz\n\n"

        with open("input.inp", "w") as inp:
            inp.writelines(input)

        return

    def spe(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 750,
        save_cubes: bool = False,
        cube_dim: int = 250,
        hirshfeld: bool = False,
        inplace: bool = False,
        remove_tdir: bool = True,
        blocks: Dict[str, Dict[str, Any]] = {},
    ):
        f"""Single point energy calculation.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 750
        save_cubes: bool, optional
            if set to True, will save a cube file containing electronic and spin densities,
            by default False.
        cube_dim: int, optional
            resolution for the cube files (default 250)
        hirshfeld: bool
            if set to true, will run the Hirshfeld population analysis. (default: False)
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        newmol : System object
            Output molecule containing the new energies.
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} SPE")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_spe",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.cube_dim = None if save_cubes is False else cube_dim
            job_info.hirshfeld = hirshfeld
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            if inplace is False:
                newmol = System(f"{mol.name}.xyz", charge=mol.charge, spin=mol.spin)
                newmol.properties = copy.copy(mol.properties)
                self.parse_output(newmol)

            else:
                self.parse_output(mol)

            process_output(mol, self.__output_suffix, "spe", mol.charge, mol.spin, save_cubes=save_cubes)

            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def opt(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 750,
        save_cubes: bool = False,
        cube_dim: int = 250,
        hirshfeld: bool = False,
        inplace: bool = False,
        remove_tdir: bool = True,
        optimization_level: Optional[str] = None,
        frequency_analysis: bool = True,
        blocks: Dict[str, Dict[str, Any]] = {},
    ):
        """Geometry optimization + frequency analysis.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 750
        save_cubes: bool, optional
            if set to True, will save a cube file containing electronic and spin densities,
            by default False.
        cube_dim: int, optional
            resolution for the cube files (default 250)
        hirshfeld: bool
            if set to true, will run the Hirshfeld population analysis. (default: False)
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        optimization_level: str
            The convergence level to be adopted during the geometry optimization (Default: NORMALOPT)
        frequency_analysis: bool
            If set to True (default) will also compute the vibration modes of the molecule and the frequencies. If the
            optimization is run in solvent, it will automatically switch to numerical frequencies.
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        newmol : System object
            Output molecule containing the new geometry and energies.
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} OPT")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_opt",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.opt = True
            job_info.freq = frequency_analysis
            job_info.cube_dim = None if save_cubes is False else cube_dim
            job_info.hirshfeld = hirshfeld
            job_info.optimization_level = optimization_level
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            if inplace is False:
                newmol = System("input.xyz", charge=mol.charge, spin=mol.spin)
                newmol.name = mol.name
                newmol.geometry.level_of_theory_geometry = self.level_of_theory
                self.parse_output(newmol)

            else:
                mol.geometry.load_xyz("input.xyz")
                mol.geometry.level_of_theory_geometry = self.level_of_theory
                self.parse_output(mol)

            process_output(mol, self.__output_suffix, "opt", mol.charge, mol.spin, save_cubes=save_cubes)

            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def opt_ts(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 750,
        save_cubes: bool = False,
        cube_dim: int = 250,
        hirshfeld: bool = False,
        inplace: bool = False,
        remove_tdir: bool = True,
        scf_convergence_level: Optional[str] = "TIGHTSCF",
        convergence_strategy: Optional[str] = "SLOWCONV",
        calculate_hessian: bool = True,
        frequency_analysis: bool = True,
        blocks: Dict[str, Dict[str, Any]] = {},
    ):
        """Transition state optimization + frequency analysis.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 750
        save_cubes: bool, optional
            if set to True, will save a cube file containing electronic and spin densities,
            by default False.
        cube_dim: int, optional
            resolution for the cube files (default 250)
        hirshfeld: bool
            if set to true, will run the Hirshfeld population analysis. (default: False)
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        scf_convergence_level: Optional[str]
            The SCF level of convergence to be adopted during the transition state optimization (Default: TIGHTSCF)
        convergence_strategy: Optional[str]
            The SCF convergence strategy to be adopted during the transition state optimization (Default: SLOWCONV)
        calculate_hessian: bool
            If set to True (default) will compute the exact Hessian before the first optimization step.
        frequency_analysis: bool
            If set to True (default) will also compute the vibration modes of the molecule and the frequencies. If the
            optimization is run in solvent, it will automatically switch to numerical frequencies.
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        newmol : System object
            Output molecule containing the new geometry and energies.
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} OPT")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_optTS",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.opt_ts = True
            job_info.freq = frequency_analysis
            job_info.cube_dim = None if save_cubes is False else cube_dim
            job_info.calc_hess = calculate_hessian
            job_info.hirshfeld = hirshfeld
            job_info.scf_convergence_level = scf_convergence_level
            job_info.scf_convergence_strategy = convergence_strategy
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            if inplace is False:
                newmol = System("input.xyz", charge=mol.charge, spin=mol.spin)
                newmol.name = mol.name
                newmol.geometry.level_of_theory_geometry = self.level_of_theory
                self.parse_output(newmol)

            else:
                mol.geometry.load_xyz("input.xyz")
                mol.geometry.level_of_theory_geometry = self.level_of_theory
                self.parse_output(mol)

            process_output(mol, self.__output_suffix, "optTS", mol.charge, mol.spin, save_cubes=save_cubes)

            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def freq(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 750,
        inplace: bool = False,
        remove_tdir: bool = True,
        raman: bool = False,
        overtones: bool = False,
        blocks: Dict[str, Dict[str, Any]] = {},
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
            memory per core, in MB, by default 750
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        raman: bool
            If set to True will compute the Raman spectrum. (default: False)
        overtones: bool
            If set to True will enable the computation of infrared overtones and combination
            bands. (default: False)
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        newmol : System object
            Output molecule containing the new energies.
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} FREQ")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_freq",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.freq = True
            job_info.raman = raman
            job_info.nearir = overtones
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            if inplace is False:
                newmol = System(f"{mol.name}.xyz", charge=mol.charge, spin=mol.spin)
                newmol.properties = copy.copy(mol.properties)
                self.parse_output(newmol)

            else:
                self.parse_output(mol)

            process_output(mol, self.__output_suffix, "freq", mol.charge, mol.spin)

            if remove_tdir:
                shutil.rmtree(tdir)

            if inplace is False:
                return newmol

    def nfreq(
        self,
        mol: System,
        ncores: int = None,
        maxcore: int = 750,
        inplace: bool = False,
        remove_tdir: bool = True,
        raman: bool = False,
        overtones: bool = False,
        blocks: Dict[str, Dict[str, Any]] = {},
    ):
        """Frequency analysis (numerical frequencies).

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 750
        inplace : bool, optional
            updates info for the input molecule instead of outputting a new molecule object,
            by default False
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        raman: bool
            If set to True will compute the Raman spectrum.
        overtones: bool
            Is set to True will enable the computation of infrared overtones and combination
            bands.
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        newmol : System object
            Output molecule containing the new energies.
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} NFREQ")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_nfreq",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.nfreq = True
            job_info.raman = raman
            job_info.nearir = overtones
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            if inplace is False:
                newmol = System(f"{mol.name}.xyz", charge=mol.charge, spin=mol.spin)
                newmol.properties = copy.copy(mol.properties)
                self.parse_output(newmol)

            else:
                self.parse_output(mol)

            process_output(mol, self.__output_suffix, "numfreq", mol.charge, mol.spin)
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
        maxcore: int = 750,
        remove_tdir: bool = True,
        blocks: Dict[str, Dict[str, Any]] = {},
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
            memory per core, in MB, by default 750
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        scan_list : Ensemble object
            Output Ensemble containing the scan frames.
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} SCAN")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_scan",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.scan = scan
            job_info.constraints = constraints
            job_info.invert_constraints = invertconstraints
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            xyz_list = [
                xyz
                for xyz in os.listdir(".")
                if os.path.splitext(xyz)[1] == ".xyz"
                and os.path.splitext(xyz)[0][:5] == "input"
                and xyz != "input.xyz"
                and xyz != "input_trj.xyz"
            ]

            mol_list = []

            # ---> evaluate if this section should/could be included in parse_output
            energies = []
            with open("output.out", "r") as f:
                read_energies = False
                for line in f:
                    if "The Calculated Surface using the SCF energy" in line:
                        read_energies = True
                        continue
                    if read_energies:
                        if len(line.split()) == 2:
                            energies.append(float(line.split()[-1]))
                        else:
                            break
            # <---

            for xyz in xyz_list:
                index = xyz.split(".")[1]
                shutil.move(f"input.{index}.xyz", f"{mol.name}.{index}.xyz")
                system = System(f"{mol.name}.{index}.xyz", charge=mol.charge, spin=mol.spin)
                system.properties.set_electronic_energy(energies.pop(0), self)
                mol_list.append(system)

            ensemble = Ensemble(mol_list)

            process_output(mol, self.__output_suffix, "scan", mol.charge, mol.spin)
            if remove_tdir:
                shutil.rmtree(tdir)

            return ensemble

    def scan_ts(
        self,
        mol: System,
        scan: str = None,
        fullscan: bool = False,
        constraints: str = None,
        invertconstraints: bool = False,
        frequency_analysis: bool = True,
        scf_convergence_level: Optional[str] = "TIGHTSCF",
        convergence_strategy: Optional[str] = "SLOWCONV",
        ncores: int = None,
        maxcore: int = 750,
        inplace: bool = False,
        remove_tdir: bool = True,
        blocks: Dict[str, Dict[str, Any]] = {},
    ):
        """
        Relaxed surface scan based transition state search.

        Parameters
        ----------
        mol : System object
            input molecule to use in the calculation
        scan : str
            string for the scan section in the %geom block
        fullscan: bool
            If set to True will not stop the scan when the transition state has been located (default: False)
        constraints : str
            string for the constraints section in the %geom block
        invertconstraints : bool, optional
            if True, treats the constraints block as the only coordinate NOT to constrain
        frequency_analysis: bool
            If set to True (default) will also compute the vibration modes of the molecule and the frequencies. If the
            optimization is run in solvent, it will automatically switch to numerical frequencies.
        scf_convergence_level: Optional[str]
            The SCF level of convergence to be adopted during the transition state optimization (Default: TIGHTSCF)
        convergence_strategy: Optional[str]
            The SCF convergence strategy to be adopted during the transition state optimization (Default: SLOWCONV)
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 750
        inplace : bool
            If set to True will update the given `mol` System with the optimized transition state structure. If set to
            False (default) will return a new system object with the optimized transition state structure together with
            an Ensamble object encoding the explored scan steps.
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction

        Returns
        -------
        System
            The optimized transition state structure. (only if inplace is False)
        Ensemble
            Output Ensemble containing the scan frames. (only if inplace is False)
        """

        logger.info(f"{mol.name}, charge {mol.charge} spin {mol.spin} - {self.method} SCAN")

        tdir = mkdtemp(
            prefix=mol.name + "_",
            suffix=f"_{self.__output_suffix}_scanTS",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):
            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if mol.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.scan_ts = scan
            job_info.freq = frequency_analysis
            job_info.fullscan = fullscan
            job_info.constraints = constraints
            job_info.invert_constraints = invertconstraints
            job_info.scf_convergence_level = scf_convergence_level
            job_info.scf_convergence_strategy = convergence_strategy
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=mol, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            if inplace is True:
                mol.geometry.load_xyz("input.xyz")
                mol.geometry.level_of_theory_geometry = self.level_of_theory
                self.parse_output(mol)

                process_output(mol, self.__output_suffix, "scanTS", mol.charge, mol.spin)
                if remove_tdir:
                    shutil.rmtree(tdir)

            else:
                xyz_list = [
                    xyz
                    for xyz in os.listdir(".")
                    if xyz.endswith(".xyz")
                    and os.path.splitext(xyz)[0][:5] == "input"
                    and xyz != "input.xyz"
                    and xyz != "input_trj.xyz"
                ]

                energies = []
                with open("output.out", "r") as f:
                    read_energies = False
                    for line in f:
                        if "The Calculated Surface using the SCF energy" in line:
                            read_energies = True
                            continue
                        if read_energies:
                            if len(line.split()) == 2:
                                energies.append(float(line.split()[-1]))
                            else:
                                break

                mol_list = []
                for xyz in xyz_list:
                    if xyz.endswith(".refined.xyz"):
                        continue

                    index = xyz.split(".")[1]
                    shutil.copy(f"input.{index}.xyz", f"{mol.name}.{index}.xyz")
                    system = System(f"{mol.name}.{index}.xyz", charge=mol.charge, spin=mol.spin)
                    system.properties.set_electronic_energy(energies[int(index) - 1], self)
                    mol_list.append(system)

                ensemble = Ensemble(mol_list)

                shutil.move("input.xyz", f"{mol.name}_TS.xyz")
                newmol = System(f"{mol.name}_TS.xyz", charge=mol.charge, spin=mol.spin)
                self.parse_output(newmol)

                process_output(mol, self.__output_suffix, "scanTS", mol.charge, mol.spin)
                if remove_tdir:
                    shutil.rmtree(tdir)

                return newmol, ensemble

    def neb_ci(
        self,
        reactant: System,
        product: System,
        nimages: Optional[int] = None,
        preoptimize: bool = False,
        ncores: int = None,
        maxcore: int = 750,
        remove_tdir: bool = True,
        blocks: Dict[str, Dict[str, Any]] = {},
    ) -> Ensemble:
        """
        Run a climbing image nudged elastic band calculation (NEB-CI) and output the ensamble encoding the optimized
        minimum energy path trajectory.

        Arguments
        ---------
        reactant: System
            The starting structure to be used in the calculation
        product: System
            The final structure to be used in the calculation
        nimages: int
            The number of images (without the fixed endpoints) to be used in the calcluation (default: 8)
        preoptimize: bool
            If set to True, will run a preoptimization in internal coordinates of the reactant and product structures.
        ncores : int, optional
            number of cores, by default all available cores
        maxcore : int, optional
            memory per core, in MB, by default 750
        remove_tdir : bool, optional
            temporary work directory will be removed, by default True
        blocks : Dict[str, Dict[str, Any]]
            The dictionary of dictionaries encoding a series of custom blocks defined by the user. If set to a non-empty
            value, will overvrite the block option eventually set on the `OrcaInput` class construction
        
        Raises
        ------
        RuntimeError
            Exception raised if the given `System` objects are not compatible with a NEB-CI calculation (e.g. same name 
            or different charge or spin multiplicity)

        Returns
        -------
        Ensamble
            The ensamble object containing the structures along the minimum energy path.        
        """

        logger.info(f"Running a NEB-CI calculation - {self.method}")
        logger.info(f"Reactant: {reactant.name}, charge {reactant.charge} spin {reactant.spin}")
        logger.info(f"Product:  {product.name}, charge {product.charge} spin {product.spin}")

        if reactant.name == product.name:
            logger.error("NEB-CI required with reactant and product with the same name")
            raise RuntimeError("Reactant and product must have different names")

        if reactant.spin != product.spin:
            logger.error("NEB-CI required with reactant and product having different spin multiplicities.")
            raise RuntimeError("Reactant and product must have the same spin multiplicity")

        if reactant.charge != product.charge:
            logger.error("NEB-CI required with reactant and product having different charge.")
            raise RuntimeError("Reactant and product must have the same charge")

        tdir = mkdtemp(
            prefix=reactant.name + "_" + product.name + "_",
            suffix=f"_{self.__output_suffix}_NEB-CI",
            dir=os.getcwd(),
        )

        with sh.pushd(tdir):

            product.geometry.write_xyz(f"{product.name}.xyz")

            job_info = OrcaJobInfo()
            job_info.ncores = ncores
            job_info.maxcore = maxcore
            job_info.is_singlet = True if reactant.spin == 1 else False
            job_info.solvent = self.solvent
            job_info.neb_ci = True
            job_info.neb_product = f"{product.name}.xyz"
            job_info.neb_images = nimages
            job_info.neb_preopt = preoptimize
            
            job_info.user_blocks = blocks if blocks != {} else self.blocks

            self.write_input(mol=reactant, job_info=job_info)

            cmd = f"{self.__ORCADIR}/orca input.inp > output.out '{cfg.MPI_FLAGS}'"
            logger.debug(f"Running Orca with command: {cmd}")
            os.system(cmd)

            MEP_systems = split_multixyz(reactant, "input_MEP_trj.xyz", suffix="MEP", engine=self)
            MEP_ensemble = Ensemble(MEP_systems) 

            if remove_tdir:
                shutil.rmtree(tdir)

            return MEP_ensemble
            

    def parse_output(self, mol: System) -> None:
        """
        The function will parse an ORCA output file automatically looking for all the relevant
        numerical properties derived form a calculation. All the properties of the given molecule
        will be set or updated.

        Parameters
        ----------
        mol: System
            The System to which the properties must be written to.

        Raises
        ------
        RuntimeError
            Exception raised if the given path to the output file is not valid.
        """

        if not os.path.isfile("output.out"):
            raise RuntimeError("Cannot parse output, the `output.out` file does not exist.")

        normal_termination = False
        with open("output.out", "r") as outfile:
            for line in outfile:
                if "****ORCA TERMINATED NORMALLY****" in line:
                    normal_termination = True
                    break

        if normal_termination is False:
            logger.error("Error occurred during orca calculation.")
            raise RuntimeError("Error occurred during orca calculation")

        # Parse the final single point energy and the vibronic energy
        # -----------------------------------------------------------------------------------
        with open("output.out", "r") as out:
            for line in out:
                if "FINAL SINGLE POINT ENERGY" in line:
                    electronic_energy = float(line.split()[-1])
                    mol.properties.set_electronic_energy(electronic_energy, self)
                if "G-E(el)" in line:
                    vibronic_energy = float(line.split()[-4])
                    mol.properties.set_vibronic_energy(vibronic_energy, self)
                if "Final Gibbs free energy" in line:
                    gibbs_free_energy = float(line.split()[-2])
                    mol.properties.set_gibbs_free_energy(gibbs_free_energy, self, self)

        # Parse the Mulliken atomic charges and spin populations
        # -----------------------------------------------------------------------------------
        counter = 0
        mulliken_charges, mulliken_spins = [], []
        spin_available = False
        with open("output.out", "r") as file:
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
                    # Check if the section contains also the "SPIN" column (either "SPIN POPULATIONS" or "SPIN DENSITIES")
                    if "SPIN" in line:
                        spin_available = True

                    _ = file.readline()  # Skip the table line

                    # Iterate over the whole section reading line by line
                    while True:
                        buffer = file.readline()
                        if "Sum of atomic charges" in buffer:
                            break
                        else:
                            data = buffer.replace(":", "").split()
                            mulliken_charges.append(float(data[2]))

                            if spin_available:
                                mulliken_spins.append(float(data[3]))
                            else:
                                mulliken_spins.append(0.0)
                else:
                    continue

                # If break has been called after mulliken has been modified the section end
                # has been reached, as such, break also from the reading operation
                if mulliken_charges != []:
                    break

        if mulliken_charges != []:
            mol.properties.set_mulliken_charges(mulliken_charges, self)
            mol.properties.set_mulliken_spin_populations(mulliken_spins, self)

        # Parse the Hirshfeld atomic charges and spin populations
        # -----------------------------------------------------------------------------------
        hirshfeld_charges, hirshfeld_spins = [], []
        with open("output.out", "r") as file:
            for line in file:
                # Read the file until the HIRSHFELD ANALYSIS title is found
                if "HIRSHFELD ANALYSIS" in line:
                    # Discard the following 6 lines to skip formatting and total integrated
                    # densities
                    for i in range(6):
                        _ = file.readline()

                    # Read the whole hirshfeld section until a empty line is found
                    while True:
                        # Read the next line
                        buffer = file.readline()

                        # If the line is empty then break else parse the line
                        if buffer == "\n":
                            break

                        else:
                            data = buffer.split()
                            hirshfeld_charges.append(float(data[2]))
                            hirshfeld_spins.append(float(data[3]))

                elif hirshfeld_charges != []:
                    break

        if hirshfeld_charges != []:
            mol.properties.set_hirshfeld_charges(hirshfeld_charges, self)
            mol.properties.set_hirshfeld_spin_populations(hirshfeld_spins, self)

        # If available parse the section related to the vibrational analysis
        # -----------------------------------------------------------------------------------
        with open("output.out", "r") as file:
            vibrational_data = None

            for line in file:
                if "VIBRATIONAL FREQUENCIES" in line:
                    vibrational_data = VibrationalData()

                    # Discard the following 4 lines to skip formatting
                    for i in range(4):
                        _ = file.readline()

                    # Read the whole vibrational frequencies section
                    while True:
                        # Read the line
                        buffer = file.readline()

                        # Break if the line is empty
                        if buffer == "\n":
                            break

                        # Strip the endline character
                        buffer = buffer.rstrip("\n")

                        # Parse the frequency line and append it to the vibrational_data class
                        substring = buffer.split(":")[-1]

                        # Check if the mode is imaginary and strip the warning
                        if substring.endswith(" ***imaginary mode***"):
                            logger.warning("Imaginary mode detected in frequency analysis.")
                            substring = substring.rstrip(" ***imaginary mode***")

                        # Parse the frequency value
                        frequency = float(substring.rstrip("cm**-1"))
                        vibrational_data.frequencies.append(frequency)

                elif "NORMAL MODES" in line:
                    # Discard the following 6 lines to skip formatting
                    for i in range(6):
                        _ = file.readline()

                    block = 0
                    while True:
                        # Discard the header line of the table block
                        _ = file.readline()

                        # Read the data within the current table block
                        ncoords = 3 * mol.geometry.atomcount
                        modes_left = ncoords - 6 * block

                        # If all the blocks have been already readed, break
                        if modes_left <= 0:
                            break

                        # Compute the number of data columns in the block
                        ncols = 6 if modes_left > 6 else modes_left

                        # Read each vector line by line
                        modes_buffer = [[] for i in range(ncols)]
                        for _ in range(ncoords):
                            sline = file.readline().split()

                            for i, element in enumerate(sline[1::]):
                                modes_buffer[i].append(float(element))

                        # Add all the obtained vectors to the vibrational data class
                        for vector in modes_buffer:
                            vibrational_data.normal_modes.append(np.array(vector))

                        # Increment the block counter
                        block += 1

                elif "IR SPECTRUM" in line:
                    # Discard the following 5 lines to skip formatting
                    for i in range(5):
                        _ = file.readline()

                    while True:
                        # Read the table line by line
                        line = file.readline()

                        # Check if the end of the table has been reached
                        if line == "\n":
                            break

                        # Split the mode index field from the rest of the data
                        sline = line.split(":")

                        # Add the mode index and the transition intensity in km/mol
                        vibrational_data.ir_transitions.append((int(sline[0]), float(sline[1].split()[2])))

                elif "OVERTONES AND COMBINATION BANDS" in line:
                    # Discard the following 5 lines to skip formatting
                    for i in range(5):
                        _ = file.readline()

                    while True:
                        # Read the table line by line
                        line = file.readline()

                        # Check if the end of the table has been reached
                        if line == "\n":
                            break

                        # Split the mode index field from the rest of the data
                        sline = line.split(":")

                        mode_index = [int(x) for x in sline[0].split("+")]
                        transition_intensity = float(sline[1].split()[2])

                        # Add the modes indeces and the transition intensity in km/mol
                        vibrational_data.ir_combination_bands.append(
                            (mode_index[0], mode_index[1], transition_intensity)
                        )

                elif "RAMAN SPECTRUM" in line:
                    # Discard the following 4 lines to skip formatting
                    for i in range(4):
                        _ = file.readline()

                    while True:
                        # Read the table line by line
                        line = file.readline()

                        # Check if the end of the table has been reached
                        if line == "\n":
                            break

                        # Split the mode index field from the rest of the data
                        sline = line.split(":")
                        mode_index = int(sline[0])
                        activity = float(sline[1].split()[1])
                        depolarization = float(sline[1].split()[2])

                        # Add the mode index, activity and depolarization
                        vibrational_data.raman_transitions.append((mode_index, activity, depolarization))

            if vibrational_data is not None:
                mol.properties.set_vibrational_data(vibrational_data, self)

    @property
    def output_suffix(self) -> str:
        """
        Suffix used to compose the name of calculation output files

        Returns
        -------
        str
            The output suffix string
        """
        return self.__output_suffix


class M06(OrcaInput):
    def __init__(self):
        super().__init__(
            method="M062X",
            basis_set="def2-TZVP",
            aux_basis="def2/J",
            solvent="water",
            optionals="DEFGRID3",
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
