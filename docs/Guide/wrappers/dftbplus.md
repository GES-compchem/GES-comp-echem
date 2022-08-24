(Guide-wrappers-dftbplus)=
# `compechem.wrappers.dftbplus` submodule

The `dftbplus` submodule contains a series of classes and functions for interfacing with the [DFTB+](https://dftbplus.org/) program.

---

## `DFTBInput` class

This submodule contains the `DFTBInput` class, which is the main interface for runnign DFTB+ calculations, specifying the Hamiltonian and its parameters.

You can initialise a `DFTBInput` wrapper by assigning an instance to a variable. In this case, the default options will be used:

```python
dftb = DFTBInput()
```

The following optional arguments are available:

* `hamiltonian` (`str`, default: `"DFTB"`): Hamiltonian used in the calculation. Available options are `DFTB`, `xTB`
* `parameters` (`str`, default: `"3ob/3ob-3-1/"`): parameters to be used for the Hamiltonian
* `solver` (`str`): LAPACK eigensolver method. If not specified, uses the internal DFTB+ default option
* `dispersion` (`bool`, default: `False`): activates D3 dispersion corrections
* `parallel` (`str`, default: `"mpi"`): selects either openmpi-parallel version (`mpi`) or shared memory version (`nompi`)
* `verbose` (`bool`, default: `False`): if `True`, saves the full DFTB+ output, otherwise, only the smaller files

:::{admonition} Note
:class: warning
For Molecular Dynamics simulations, the output files can become very large (several hundred GB in some cases), make sure you have enough space on your disk if you activate the `verbose` option!
:::

Please consult the DFTB+ manual for all available options and settings.

---

## Methods

The `DFTBInput` class features a series of internal methods to carry out the calculations. The available options are:

* `spe`: single point energy calculation
* `md_nvt`: molecular dynamics simulation in the canonical (NVT) ensemble
* `simulated_annealing`: simulated annealing simulation 

These methods can be called directly from the `DFTBInput` instance, passing the `System` object on which you want to carry out the calculation:

```python
mymol = System("mymol.xyz")
dftb = DFTBInput()

mymol_opt = dftb.spe(mymol)
```

---

## `DFTBInput.spe()` parameters

The single point energy function accepts the following optional flags:

* `ncores` (`int`): number of cores for the calculation. By default takes all available cores on the machine
* `charge` (`int`): charge of the system. By default, taken from the input system
* `charge` (`int`): spin multiplicity of the system. By default, taken from the input system
* `inplace` (`bool`, default: `False`): if `True`, updates the information (energies, geometry, etc.) of the original `System`. Otherwise, returns a new `System` which can be assigned to a new variable
* `remove_tdir` (`bool`, default: `True`): if `True`, deletes the temporary work directory after a successful calculation. 

---

## `DFTBInput.md_nvt()` parameters

The NVT MD function **requires** the `steps` parameter, setting the number of MD simulation steps, as well as accepting the following optional flags:

* `timestep` (`float`, default: `1.0`): time step (in fs)
* `temperature` (`float`, default: `298.0`): temperature (in Kelvin)
* `mdrestartfreq` (`int`, default: `100`): data is saved every mdrestartfreq steps, by default 100
* `box_side` (`float`): for periodic systems, defines the length (in Å) of the box side. By default, taken from the input `System`
* `ncores` (`int`): number of cores for the calculation. By default takes all available cores on the machine
* `charge` (`int`): charge of the system. By default, taken from the input system
* `charge` (`int`): spin multiplicity of the system. By default, taken from the input system
* `inplace` (`bool`, default: `False`): if `True`, updates the information (energies, geometry, etc.) of the original `System`. Otherwise, returns a new `System` which can be assigned to a new variable
* `remove_tdir` (`bool`, default: `True`): if `True`, deletes the temporary work directory after a successful calculation. 
* `compress_traj` (`bool`, default: `True`): if True, parses the geo.end and md.out files into a single, smaller file, which is then zipped in an archive.

---

## `DFTBInput.simulated_annealing()` parameters

The Simulated Annealing accepts the following optional flags:

* `start_temp` (`float`, default: `1.0`): starting temperature in Kelvin
* `target_temp` (`float`, default: `2000.0`): maximum temperature in Kelvin
* `ramp_steps` (`int`, default: `500`): number of MD steps for the heating/cooling ramps
* `hold_steps` (`int`, default: `1000`): number of MD steps held at target_temp
* `timestep` (`float`, default: `1.0`): time step (in fs)
* `mdrestartfreq` (`int`, default: `100`): data is saved every mdrestartfreq steps, by default 100
* `box_side` (`float`): for periodic systems, defines the length (in Å) of the box side. By default, taken from the input `System`
* `ncores` (`int`): number of cores for the calculation. By default takes all available cores on the machine
* `charge` (`int`): charge of the system. By default, taken from the input system
* `charge` (`int`): spin multiplicity of the system. By default, taken from the input system
* `inplace` (`bool`, default: `False`): if `True`, updates the information (energies, geometry, etc.) of the original `System`. Otherwise, returns a new `System` which can be assigned to a new variable
* `remove_tdir` (`bool`, default: `True`): if `True`, deletes the temporary work directory after a successful calculation. 
* `compress_traj` (`bool`, default: `True`): if True, parses the geo.end and md.out files into a single, smaller file, which is then zipped in an archive.
