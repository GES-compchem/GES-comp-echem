(Guide-wrappers-orca)=
# `compechem.wrappers.orca` submodule

The `orca` submodule contains a series of classes and functions for interfacing with the [Orca](https://orcaforum.kofo.mpg.de/app.php/dlext/) program.

---

### `OrcaInput` class

This submodule contains the `OrcaInput` class, which is the main interface for runnign Orca calculations, specifying the level of theory at which calculations will be carried out, the basis set to be used, if and which implicit solvent should be used, etc.

You can initialise an `OrcaInput` wrapper by assigning an instance to a variable. In this case, the default options will be used:

```python
orca = OrcaInput()
```

The following optional arguments are available:

* `method` (`str`): Level of theory (or XC functional) for the calculation
* `basis_set` (`str`, default: `def2-TZVP`): basis set to be used in the calculation
* `aux_basis` (`str`, default: `def2/J`): auxiliary basis set to be used in the calculation
* `solvation` (`bool`, default: `False`): activates the CPCM(SMD) implicit solvation model
* `solvent` (`str`, default: `"water"`): SMD solvent
* `optionals` (`str`, default: `""`): optional keywords/flags to be included in the Orca input file

Please consult the Orca manual for all available options for methods and basis sets.

---

### Methods

The `OrcaInput` class features a series of internal methods to carry out the calculations. The available options are:

* `spe`: single point energy calculation
* `opt`: geometry optimisation + frequency analysis
* `freq`: analytical frequency analysis
* `nfreq`: numerical frequency analysis only

These methods can be called directly from the `OrcaInput` instance, passing the `System` object on which you want to carry out the calculation:

```python
mymol = System("mymol.xyz")
orca = OrcaInput()

mymol_opt = orca.opt(mymol)
```

All the methods above accept a series of optional flags:

* `ncores` (`int`): number of cores for the calculation. By default takes all available cores on the machine
* `maxcore` (`int`, default: `350`): maximum available memory (in MB), per core
* `charge` (`int`): charge of the system. By default, taken from the input system
* `spin` (`int`): spin multiplicity of the system. By default, taken from the input system
* `inplace` (`bool`, default: `False`): if `True`, updates the information (energies, geometry, etc.) of the original `System`. Otherwise, returns a new `System` which can be assigned to a new variable
* `remove_tdir` (`bool`, default: `True`): if `True`, deletes the temporary work directory after a successful calculation. 

---

### Pre-packaged methods

For some common types of calculation, some "pre-packaged" options have been set up and made available:

* `M06`: M06-2X density functional with the def2-TZVP basis set and D3ZERO dispersion corrections, in implicit SMD water solvent
* `r2SCAN`: r2SCAN-3c density functional with D3ZERO dispersion corrections, in implicit SMD water solvent
* `CCSD`: DLPNO-CCSD calculation extrapolating to the ANO basis set limit (using double-zeta and triple-zeta points) in implicit SMD water solvent

These can be called as a normal `OrcaInput` instance, directly assigning it to a variable:

```python
from compchem import M06
m06 = M06()