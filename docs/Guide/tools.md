(Guide-tools)=
# `compechem.tools` submodule

The `tools` submodule contains a series of useful functions, mostly used internally by the library. The user is not expected to interact with this submodule during normal operation, but some functions might come in handy for debugging or advanced users, and will be reported here. For a comprehensive list of all functions available in the `tools` submodule, including the ones only meant for internal use, please refer to the [API](API-tools).

The `tools` module can be imported via the following syntax:

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
