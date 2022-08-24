(Guide-wrappers-xtb)=
# `compechem.wrappers.xtb` submodule

The `xtb` submodule contains a series of classes and functions for interfacing with the [xTB](https://github.com/grimme-lab/xtb) program.

---

### `XtbInput` class

This submodule contains the `XtbInput` class, which is the main interface for runnign xTB calculations, specifying the level of theory at which calculations will be carried out, if and which implicit solvent should be used, etc.

You can initialise an `XtbInput` wrapper by assigning an instance to a variable. In this case, the default options will be used:

```python
xtb = XtbInput()
```

The following optional arguments are available:

* `method` (`str`, default: `"gfn2"`): Hamiltonian used in the calculation. Available options are `gfn2`, `gfn1`, `gfnff`
* `solvation` (`bool`, default: `True`): activates the ALPB implicit solvation model
* `solvent` (`str`, default: `"water"`): ALPB solvent
* `optionals` (`str`, default: `""`): optional keywords/flags to be passed to xTB command line call. Consult the [xTB](https://xtb-docs.readthedocs.io/en/latest/contents.html) manual for a list of all the available options

---

### Methods

The `XtbInput` class features a series of internal methods to carry out the calculations. The available options are:

* `spe`: single point energy calculation
* `opt`: geometry optimisation + frequency analysis
* `freq`: frequency analysis only

These methods can be called directly from the `XtbInput` instance, passing the `System` object on which you want to carry out the calculation:

```python
mymol = System("mymol.xyz")
xtb = XtbInput()

mymol_opt = xtb.opt(mymol)
```

All the methods above accept a series of optional flags:

* `ncores` (`int`): number of cores for the calculation. By default takes all available cores on the machine
* `charge` (`int`): charge of the system. By default, taken from the input system
* `spin` (`int`): spin multiplicity of the system. By default, taken from the input system
* `inplace` (`bool`, default: `False`): if `True`, updates the information (energies, geometry, etc.) of the original `System`. Otherwise, returns a new `System` which can be assigned to a new variable
* `remove_tdir` (`bool`, default: `True`): if `True`, deletes the temporary work directory after a successful calculation. 