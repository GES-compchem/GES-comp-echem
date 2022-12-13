import subprocess
from tempfile import NamedTemporaryFile as tmp

from os import system
from os.path import join, abspath, basename


def locate_vmd() -> str:
    """
    Locate the path to 'vmd' from the system PATH.

    Raises
    ------
    RuntimeError
        Exception raised if vmd is not found

    Returns
    -------
    str
        The path to the vmd base folder (the one containing the `bin` and `lib` subfolders)
    """
    output = subprocess.check_output("whereis vmd", shell=True).decode("utf-8").strip("\n")

    if len(output.split(": ")) == 1:
        raise RuntimeError("cannot find 'vmd' in the system path")
    else:
        return abspath(output.split(": ")[-1].rstrip("/bin/vmd"))


def render_fukui_cube(
    cubfile: str,
    isovalue: float = 0.003,
    include_negative: bool = False,
    resolution: int = 4800,
    shadows: bool = True,
    ambientocclusion: bool = True,
    dof: bool = True,
    VMD_PATH: str = None,
) -> None:
    """
    Given the path to a Fukui function cube file saves a `.png` render of the positive part
    of the function.

    Arguments
    ---------
    cubefile: str
        The path to the `.fukui.cube` file that must be rendered.
    isovalue: float
        The isovalue at which the contour must be plotted (default: 0.003).
    include_negative: bool
        If set to True, will render also the negative part of the Fukui function. (default:
        False)
    resolution: int
        The resolution of the output image (default: 4800).
    shadows: bool
        If set to True will enable the vmd shadows option
    ambientocclusion: bool
        If set to True will enable the vmd ambientocclusion option
    dof: bool
        If set to True will enable the vmd dof option
    VMD_PATH: str
        The path to the vmd folder. Is set to None (default), will automatically search vmd
        in the system PATH.
    """
    vmd_root = VMD_PATH if VMD_PATH is not None else locate_vmd()
    tachyon_path = join(vmd_root, "lib/vmd/tachyon_LINUXAMD64")

    root_name = basename(cubfile).strip(".fukui.cube")

    with tmp(mode="w+") as vmd_script:

        vmd_script.write(
            f"""
            mol addrep 0
            display projection Orthographic
            display resetview
            mol new {cubfile} type {{cube}} first 0 last -1 step 1 waitfor 1 volsets {{0 }}
            animate style Loop
            axes location Off
            mol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000
            mol color Name
            mol representation CPK 1.000000 0.300000 150.000000 12.000000
            mol selection all
            mol material Opaque
            mol addrep 0
            mol modcolor 1 0 Volume 0
            mol modstyle 1 0 Isosurface {isovalue} 0 0 0 1 1
            mol modmaterial 1 0 Translucent
            mol scaleminmax 0 1 0.000000 1.000000
            """
        )

        if include_negative:
            vmd_script.write(
                f"""
                mol addrep 0
                mol modcolor 2 0 Volume 0
                mol modstyle 2 0 Isosurface {-isovalue} 0 0 0 1 1
                mol modmaterial 1 0 Translucent
                mol scaleminmax 0 2 -1.000000 0.000000
                """
            )

        if shadows:
            vmd_script.write("display shadows on\n")

        if ambientocclusion:
            vmd_script.write("display ambientocclusion on\n")

        if dof:
            vmd_script.write("display dof on")

        vmd_script.write(
            f"""
            color Display Background white
            color Element C black
            mol modcolor 0 0 Element
            render Tachyon {root_name}.dat "{tachyon_path}" -fullshade -aasamples 12 %s -format TARGA -res {resolution} {resolution} -o %s.png
            exit
            """
        )

        vmd_script.seek(0)
        system(f"vmd -dispdev text -e {vmd_script.name}")
