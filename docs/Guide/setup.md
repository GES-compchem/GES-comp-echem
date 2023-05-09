(Guide-setup)=
# Setting up calculations

After creating a `System`, you need to set up the calculation to which you wish to subject it. This is done with the `compechem.engines` and `compechem.wrappers` submodules.

Both `engines` and `wrappers` submodules contains a series of program-specific classes and functions for interfacing with external code and carrying out calculations on `System` objects. The distinction between the two and the philosophy behind them is as follows:

* an `engine` carries out calculations on a `System` for computing properties. For example, geometry optimizations, single point energies, frequencies, etc. 
* a `wrapper` carries out calculations on a `System` for obtaining other, "processed" `System` objects. For example, conformer/tautomer searches, building of solvation boxes, etc.

The programs currently implemented as `engines` are:

* [xTB](https://github.com/grimme-lab/xtb)
* [Orca](https://orcaforum.kofo.mpg.de/index.php?sid=3c6c78cae3dd0cfffa26a293953422e3)
* [DFTB+](https://dftbplus.org/)
* [NAMD](http://www.ks.uiuc.edu/Research/namd/) (coming soon!)

The programs currently implemented as `wrappers` are:

* [CREST](https://crest-lab.github.io/crest-docs/)
* [PackMol](https://m3g.github.io/packmol/)

:::{admonition} Note
:class: warning
To function with the library, the external programs need to be already installed and available to the system from command line!
:::

These general-purpose `engines` are implemented as `<Program>Input` classes in the corresponding submodules. To initialise an `engine`, you need to import it from its submodule and then create an instance of it:

```python
from compechem.engines.xtb import XtbInput
from compechem.engines.dftbplus import DFTBInput
from compechem.engines.orca import OrcaInput

xtb = XtbInput()
dftb = DFTBInput()
orca = OrcaInput()
```

If you do not specify anything, some default options are chosen automatically for level of theory, basis set, solvation, etc. Please refer to the [API](API-engines) section for a complete list of options and default values.

As an example, let us set up a calculation with Orca, using the B3LYP functional, with the def2-TZVP basis set and def2/J auxiliary basis set, using the SMD implicit solvation model for water, and including Grimme's D3BJ dispersion corrections:

```python
from compechem.engines.orca import OrcaInput

b3lyp = OrcaInput(
    method = "B3LYP",
    basis_set = "def2-TZVP",
    aux_basis = "def2/J",
    solvent = "water",
    optionals = "D3BJ",
)
```

For `wrappers`, you just need to import the corresponding `wrapper` submodule or any of the specific functions you wish to use:

```python
from compechem.wrappers.crest import tautomer_search

mol = System("water.xyz")
tautomers_list = tautomer_search(mol)
```