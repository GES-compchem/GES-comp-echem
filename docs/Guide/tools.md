(Guide-tools)=
# Useful tools

The `compechem.tools` submodule contains a series of useful functions, mostly used internally by the library. The user is not expected to interact with this submodule during normal operation, but some functions might come in handy for debugging or for advanced users, and will be reported here. For a comprehensive list of all functions available in the `compechem.tools` submodule, including the ones only meant for internal use, please refer to the [API](API-tools).

The `compechem.tools` module can be imported via the following syntax:

```python
from compechem import tools
```

while individual functions can be imported as:

```python
from compechem.tools import split_multixyz
```

---

## `cyclization_check()`

Checks if a cyclization has occurred (e.g., during a geometry optimization), or if a ring opening has occurred. 

The function accepts two arguments, `start_file` and `end_file`, both strings, containing the path to the two geometries to be compared.

The function uses the [CREST](https://github.com/grimme-lab/crest) program to generate the topology for each structure, and counts the number of rings. If a mismatch in the number of rings is observed, returns `True`.

---

## `dissociation_check()`

Checks if a dissociation has occurred (e.g., during a geometry optimization). 

The function looks for a `.mol` file in the current directory, and generates the corresponding [SMILES](https://it.wikipedia.org/wiki/SMILES) string. If a dot (".") is present in the SMILES string, it means there are two separate molecules in the structure, and therefore a dissociation has occurred.

---

## `split_multixyz()`

Splits a .xyz file containing multiple structures into individual `System` objects. For example, used to generate each separate conformer from a "global" .xyz file containing all the structures out of a CREST conformer search. 

The function requires the `System` object which will be used as "template" for the output, and the .xyz file in which the structures are collected. It returns a list of `System` objects, corresponding to the individual structures.

```python
mymol = System("all_conformers.xyz")
multi_xyz_list = tools.split_multixyz(mol=mymol, file="all_conformers.xyz", suffix="conf")
```

---

## `maxdist()`

Calculates the maximum distance (in Å) between any 2 atoms in the structure of the provided `System` object.

---

## `reorder_energies()`

Reorders a `System` list at a specific level of theory, provided the following arguments:

* `molecule_list` (`list(System)`): list containing the System objects to be reordered (e.g., generated by a CREST routine)
* `ncores` (`int`): number of cores used in the calculations, by default all available cores
* `maxcore` (`int`, default: `350`): memory per core, in MB
* `method_opt` (any Input type from a `wrapper`): level of theory for the geometry optimization. By default GFN2-xTB
* `method_el` (any Input type from a `wrapper`): level of theory for the electronic component of the total energy. By default M06-2X/D3ZERO
* `method_vib` (any Input type from a `wrapper`): level of theory for the vibronic component of the total energy. By default GFN2-xTB

and returns the same list, but with a new ordering, given by the total energies recalculated at the new level of theory.

## VMD based tools

The `copechem.tools` module also provides a `vmdtools` submodule containing simple VMD-based functions to render graphical representations of molecular properties such as:

* Volumetric and condensed Fukui functions
* Spin densities

The following conventions have been adopted in the representation:

* Positive/high values of functions/densities related to the electronic density are represented in red while negative/low values are represented in blue.
* Positive/high values of functions/densities related to the spin density are represented in orange while negative/low values are represented in violet.

Please refer to the API for details about the call to the various render functions.


## `Mogli` based tools

The `copechem.tools` module also provides a `moglitools` submodule containing simple mogli-based molecular viewer. Please refer to the API for details about the call to the various render functions.