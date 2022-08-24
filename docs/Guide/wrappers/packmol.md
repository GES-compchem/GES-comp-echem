(Guide-wrappers-packmol)=
# `compechem.wrappers.packmol` submodule

The `packmol` submodule contains a series of functions for interfacing with the [packmol](http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml) program, and can be imported with the following syntax:

```python
from compechem import packmol
```

---

### `packmol.packmol_cube()`

the `packmol_cube` function generates a cubic solvation box of **solvent** molecules around a **solute** object. It returns a `System` object corresponding to the (periodic) solvated box.

Two out of the following three parameters need to be passed to the function to correctly generate the box:

* `nsolv` (`int`): number of solvent molecules to add in the box
* `target_dens` (`float`): target density for the solvated box, in g/L
* `cube_side` (`int`): length of the box side, in Å

The third parameter will be calculated from the other two.

For example, to generate a solvated box with 100 solvent molecules and a density of 997 g/L:

```python
solute = System("solute.xyz")
solvent = System("solvent.xyz")

solvated_box = packmol.packmol_cube(solute, solvent, nsolv=100, target_dens=997)
```
