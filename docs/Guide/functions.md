(Guide-functions)=
# `compechem.functions` submodule

The `functions` submodule contains a series of methods to "manually" calculate some physical properties of the system under study. The user is expected to provide the exact intended states for the system(s) under study, as the functions will simply take the inputs and generate the output without carrying out any checks. For example, if the user calculates the pKa of a molecule in an unoptimised geometry, the resulting pKa will be that of the unoptimised geometry.

The `functions` module can be imported via the following syntax:

```python
from compechem import functions
```

while individual methods can be imported as:

```python
from compechem.functions import calculate_pka
```

---

## `calculate_pka()`

Calculates the pKa of a system, provided the following arguments:

* `protonated` (`System`): molecule in its protonated form
* `deprotonated` (`System`): molecule in its deprotonated form
* `method_el` (`System`): level of theory for the electronic component of the total energy (must be present in the `System.energies` dictionary)
* `method_vib` (`System`): level of theory for the vibronic component of the total energy (must be present in the `System.energies` dictionary). If not specified, defaults to the same level of theory as the electronic energy

and returns the pKa of the molecule considering the provided states.

---

## `calculate_potential()`

Calculates the reduction potential of a molecule, provided the following arguments

* `oxidised` (`System`): molecule in its oxidised state
* `reduced` (`System`): molecule in its reduced state
* `method_el` (`System`): level of theory for the electronic component of the total energy (must be present in the `System.energies` dictionary)
* `method_vib` (`System`): level of theory for the vibronic component of the total energy (must be present in the `System.energies` dictionary). If not specified, defaults to the same level of theory as the electronic energy
* `pH` (`float`, default: `7.0`): pH at which the reduction potential is calculated

and returns the reduction potential of the molecule considering the provided states at the provided pH, including eventual PCET mechanisms


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