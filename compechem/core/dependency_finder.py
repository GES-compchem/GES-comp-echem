import subprocess

from os.path import abspath
from io import BytesIO


def locate_executable(name: str) -> str:
    """
    Locate in the system PATH a given executable. If the program is found the path is
    returned, else an exception is raised.

    Arguments
    ---------
    name: str
        The name of the program to locate.

    Raises
    ------
    RuntimeError
        Exception raised if the program is not found in the system path.

    Returns
    -------
    str
        The path to the requested program.
    """
    output = subprocess.check_output(f"whereis {name}", shell=True).decode("utf-8")
    output = output.strip("\n")

    if len(output.split(": ")) == 1:
        raise RuntimeError(f"cannot find '{name}' in the system path")
    else:
        return str(abspath(output.split(": ")[-1]))


def locate_vmd() -> str:
    """
    Locate the path to the 'vmd' folder from the system PATH.

    Returns
    -------
    str
        The path to the vmd base folder (the one containing the `bin` and `lib` subfolders)
    """
    path = locate_executable("vmd")
    return path.rstrip("/bin/vmd")


def locate_orca(version: str = None, get_folder: bool = False) -> str:
    """
    Locate the path to the 'orca' executable from the system PATH. If the executable is
    located check that the correct version of OpenMPI is exported (explicit reference to
    the static builds). If specified, check that the correct version of orca is loaded.

    Arguments
    ---------
    version: str
        The string defining the desired version of orca. If set to None (default) all version
        of orca are accepted.
    get_folder: bool
        If set to True will return the path of the folder containing the orca executable instead
        of the default path to the executable itself. Equivalent to applying an `rstrip('/orca')`
        to the executable path
        
    Returns
    -------
    str
        The path to the orca executable file.
    """
    path = locate_executable("orca")

    # Check if the available version of orca matches the requirements
    orca_version = None
    output = subprocess.run(["orca", "--version"], capture_output=True, text=True).stdout
    for line in output.split("\n"):

        if "Program Version" in line:
            orca_version: str = line.split()[2]

        elif orca_version is not None:
            break

    if orca_version is None:
        raise RuntimeError("Failed to read the version of the orca software.")

    elif version is not None and orca_version != version:
        raise RuntimeError(
            f"The required orca version is not available. Version {orca_version} found instead."
        )

    # Check if a MPI version is available in the system PATH
    _ = locate_executable("mpirun")

    # Check if the version of the loaded OpenMPI
    openmpi_version = None
    output = subprocess.run(["mpirun", "--version"], capture_output=True, text=True).stdout

    for line in output.split("\n"):

        if "(Open MPI)" in line:
            openmpi_version: str = line.split()[3]

        elif openmpi_version is not None:
            break

    if openmpi_version is None:
        raise RuntimeError("OpenMPI is either not available or the version cannot be found")

    # Check if the available version meets the requirements
    openmpi_required ={
        "5.0.*": ["4.1.1"],
        "4.2.*": ["3.1.4"],
        "4.1.*": ["3.1.3", "2.1.5"]
    }

    key = ".".join(orca_version.split(".")[0:-1] + ["*"])
    if openmpi_version not in openmpi_required[key]:
        msg = " or ".join(openmpi_required[key])
        raise RuntimeError(
            f"orca {orca_version} retuires OpenMPI {msg}. OpenMPI {openmpi_version} found instead."
        )

    return path.rsplit("/orca") if get_folder is True else path

