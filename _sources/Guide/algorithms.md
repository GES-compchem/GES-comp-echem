(Guide-algorithms)=
# `compechem.algorithms` submodule

The `algorithms` submodule provides a series of high-level algorithms to automatically calculate some physical observables, considering the most appropriate states for the input systems. The user is only required to provide a generic starting structure for the system of interest, and the algorithm will make sure to generate the correct "state", optimising components such as geometry, conformer/tautomer, protonation state at the pH in question, etc.

The `algorithms` module can be imported via the following syntax:

```python
from compechem import algorithms
```

while individual methods can be imported as:

```python
from compechem.algorithms import one_electron_oxidation_potentials
```

---

## `one_electron_oxidation_potentials()`

Calculates the 1-electron oxidation potential of a system, provided the following arguments:

* `system` (`System`): starting system, in the reduced state
* `method` (any Input type from a `wrapper`): level of theory at which the potential is calculated.
* `ncores` (`int`): number of cores used in the calculations, by default all available cores
* `maxcore` (`int`, default: `350`): memory per core, in MB
* `conformer_search` (`bool`, default: `True`): if True, also considers the lowest-energy conformers at each stage
* `tautomer_search` (`bool`, default: `True`): if True, also considers the lowest-energy tautomer at each stage
* `pH_step` (`float`, default: `1.0`): calculates the redox potential at every pH_step units of pH

and returns an iterator which yields the 1-electron oxidation potential at each `pH_step` unit of pH.

The algorithm (optionally) determines the lowest-energy conformer and tautomer for every species under consideration. It also determines the correct protonation states, for both reduced and oxidised forms, at every pH value, optimising the geometry and calculating the energies at the requested level of theory.