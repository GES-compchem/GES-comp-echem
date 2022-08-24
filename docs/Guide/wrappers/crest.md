(Guide-wrappers-crest)=
# `compechem.wrappers.crest` submodule

The `crest` submodule contains a series of functions for interfacing with the [CREST](https://github.com/grimme-lab/crest) program, and can be imported with the following syntax:

```python
from compechem import crest
```

All functions share some common optional arguments:

* `ncores` (`int`): number of cores for the calculation. By default takes all available cores on the machine
* `remove_tdir` (`bool`, default: `True`): if `True`, deletes the temporary work directory after a successful calculation. 
* `optionals` (`str`, default: `""`): optional keywords/flags to be passed to CREST command line call. Consult the [CREST](https://xtb-docs.readthedocs.io/en/latest/crest.html) manual for a list of all the available options

---

### `crest.tautomer_search()`

the `tautomer_search` function carries out a tautomer search, taking a `System` object as input, and returning a `list` of `System` objects corresponding to the lowest energy tautomers, ordered by increasing energy:

```python
mol = System("mymol.xyz")

# all tautomers
tautomers = crest.tautomer_search(mol)

# lowest-energy tautomer
best_tautomer = tautomers[0]
```

---

### `crest.conformer_search()`

the `conformer_search` function carries out a conformer search, taking a `System` object as input, and returning a `list` of `System` objects corresponding to the lowest energy conformers, ordered by increasing energy:

```python
mol = System("mymol.xyz")

# all conformers
conformers = crest.conformer_search(mol)

# lowest-energy conformer
best_conformer = conformers[0]
```

---

### `crest.deprotonate()`

the `deprotonate` function carries out a deprotomer search, taking a `System` object as input, and returning a `list` of `System` objects corresponding to the lowest energy deprotomers, ordered by increasing energy:

```python
mol = System("mymol.xyz")

# all deprotomers
deprotomers = crest.deprotonate(mol)

# lowest-energy deprotomer
best_deprotomer = deprotomers[0]
```

---

### `crest.protonate()`

the `protonate` function carries out a protomer search, taking a `System` object as input, and returning a `list` of `System` objects corresponding to the lowest energy protomers, ordered by increasing energy:

```python
mol = System("mymol.xyz")

# all protomers
deprotomers = crest.protonate(mol)

# lowest-energy protomer
best_deprotomer = protomers[0]
```

---

### `crest.qcg_grow()`

the `qcg_grow` function carries out a [Quantum Cluster Growth](https://xtb-docs.readthedocs.io/en/latest/crestqcg.html) calculation, taking two `System` objects as input, namely a **solute** and a **solvent**, and returning another `System` corresponding to the explicitly solvated "cluster" obtained by surrounding the **solute** with **solvent** molecules:

```python
solute = System("solute.xyz")
solvent = System("solvent.xyz")

cluster = crest.qcg_grow(solute, solvent)
```

The `qcg_grow` also takes the following optional arguments:

* `charge` (`int`): charge of the system. By default, taken from the solute system
* `spin` (`int`): spin multiplicity of the system. By default, taken from the input system
* `method` (`str`, default: `"gfn2"`): Hamiltonian used in the calculation. Available options are `gfn2`, `gfn1`, `gfnff`
* `nsolv` (`int`): number of solvent molecules to put in the cluster. By default, the program keeps adding solvent molecules until energy is converged with respect to the number of solvent molecules
---

### `crest.qcg_ensemble()`

the `qcg_ensemble` function carries out a [Quantum Cluster Growth](https://xtb-docs.readthedocs.io/en/latest/crestqcg.html) calculation + an ensemble search, taking two `System` objects as input, namely a **solute** and a **solvent**, and returning an `Ensemble` corresponding to the explicitly solvated "clusters" obtained by surrounding the **solute** with **solvent** molecules, in order of ascending energy:

```python
solute = System("solute.xyz")
solvent = System("solvent.xyz")

cluster = crest.qcg_ensemble(solute, solvent)
```

The `qcg_ensemble` also takes the following optional arguments:

* `charge` (`int`): charge of the system. By default, taken from the solute system
* `spin` (`int`): spin multiplicity of the system. By default, taken from the input system
* `method` (`str`, default: `"gfn2"`): Hamiltonian used in the calculation. Available options are `gfn2`, `gfn1`, `gfnff`
* `enslvl` (`str`, default: `"gfn2"`): Hamiltonian used in the ensemble generation/evaluation. Available options are `gfn2`, `gfn1`, `gfnff`
* `ensemble_choice` (`str`, default: `"full_ensemble"`): file containing the chosen ensemble after generation. Available options are `full_ensemble`, `final_ensemble`, `crest_best`
* `nsolv` (`int`): number of solvent molecules to put in the cluster. By default, the program keeps adding solvent molecules until energy is converged with respect to the number of solvent molecules

