(Guide-functions)=
# Analysing results

The `compechem.functions` submodule contains a series of methods to "manually" calculate some physical properties of the system under study. The user is expected to provide the exact intended states for the system(s) under study, as the functions will simply take the inputs and generate the output without carrying out any checks. For example, if the user calculates the pKa of a molecule in an unoptimised geometry, the resulting pKa will be that of the unoptimised geometry.

The `compechem.functions` functions can be imported via the following syntax:

```python
from compechem.functions import calculate_pka
```

---

## pKa

The `compechem.functions.calculate_pka` function calculates the pKa of a molecule $HA$ considering a reaction of the type:

$$
HA \rightarrow H^{+} + A^{-}
$$

provided the following arguments:

* `protonated` (`System`): molecule in its protonated form
* `deprotonated` (`System`): molecule in its deprotonated form
* `method_el` (`System`): level of theory for the electronic component of the total energy (must be present in the `System.energies` dictionary)
* `method_vib` (`System`): level of theory for the vibronic component of the total energy (must be present in the `System.energies` dictionary). If not specified, defaults to the same level of theory as the electronic energy

The function returns the pKa of the molecule considering the provided states, calculated as:

$$
pK_{a} = \frac{G_{A^{-}} + G_{H^{+}} - G_{HA}}{2.303 \cdot RT}
$$

Where $G_{A^{-}}$ and $G_{HA}$ are calculated summing the electronic + vibronic energies at the selected level of theory, $G_{H^{+}} = -270.29 kcal/mol$, $R = 1.987 \cdot 10^{-3} kcal/(mol \cdot K)$ and $T = 298.15 K$

---

## 1-el redox potential

Calculates the one-electron reduction potential of a molecule $MH_{n}$, considering a generic reaction of the type:

$$
MH_{n} \rightarrow M^{\cdot (n-1)-} + n\cdot H^{+} + e^{-}
$$

provided the following arguments:

* `oxidised` (`System`): molecule in its oxidised state
* `reduced` (`System`): molecule in its reduced state
* `method_el` (`System`): level of theory for the electronic component of the total energy (must be present in the `System.energies` dictionary)
* `method_vib` (`System`): level of theory for the vibronic component of the total energy (must be present in the `System.energies` dictionary). If not specified, defaults to the same level of theory as the electronic energy
* `pH` (`float`, default: `7.0`): pH at which the reduction potential is calculated

and returns the reduction potential of the molecule considering the provided states at the provided pH, including eventual PCET mechanisms, calculated as:

$$
E°_{MH_{n}/M^{\cdot(n-1)-}} = - \frac{G_{M^{\cdot(n-1)-}} - (G_{MH_{n}} + n \cdot G_{H^{+}}) }{F} - E°_{SHE} - n \cdot 0.059 \cdot pH 
$$

where $G_{M^{\cdot(n-1)-}}$ and $G_{MH_{n}}$ are calculated summing the electronic + vibronic energies at the selected level of theory, $G_{H^{+}} = -270.29 kcal/mol$, $F = 23.061 kcal/volt–gram-equivalent$, and $E_{SHE} = 4.28 V$.

## Fukui functions
The `compechem.functions.calculate_fukui` function calculates the Fukui functions $f^+(r)$, $f^-(r)$ and $f^0(r)$ associated with a given molecular geometry. The Fukui functions are computed according to the definitions:

$$
f^+(r) = \rho_{N+1}(r) - \rho_{N}(r)
$$

$$
f^-(r) = \rho_{N}(r) - \rho_{N-1}(r)
$$

$$
f^0(r) = \frac{1}{2} \left[\rho_{N+1}(r) - \rho_{N-1}(r)\right]
$$

Where, given a molecule with $N$ electrons, $\rho_{N}(r)$ represents its electronic density while $\rho_{N\pm1}(r)$ represents the electronic density of the molecule, in the same nuclear configuration, when one electron is either added ($+1$) or removed ($-1$).

The Fukui functions are both computed as volumetric quantities and saved in [Gaussian Cube format](http://paulbourke.net/dataformats/cube/) in the `output_density` folder and as localized values returned in the form of a dictionary. The localized Fukui functions are computed by applying the $f^+$, $f^-$ and $f^0$ definitions replacing the charge density with the Mulliken charges.

The function can be called with the following minimal arguments:
* `molecule` (`System`): The molecular structure to be used in the computation
* `orca` (`OrcaInput`): The ORCA wrapper defining the level of theory to be used in the calculation.

The functions assumes that the molecule supports only singlet and doublet states and switches the spin multeplicity according to the number of electrons. If different spin states needs to be considered the `spins_states` option can be used to provide the spin multeplicity values as a list.

An example code snippet is provided in what follows:

```python
molecule = System(xyzfile)
wrapper = OrcaInput(method="B3LYP", basis_set="6-311++G**")

localized_fukui = calculate_fukui(molecule, wrapper)

print("Localized Fukui functions:")
for level_of_theory, data in localized_fukui.items():
    print(level_of_theory)
    for key, mylist in data.items():
        print(f" -> {key}: {mylist}")
```

That for a water molecule returns the following result:
```
Localized Fukui functions:
B3LYP|6-311++G**
 -> f+: [0.39742500000000003, -0.699236, -0.698188]
 -> f-: [-0.740002, -0.129992, -0.130005]
 -> f0: [-0.17128849999999998, -0.414614, -0.4140965]
```

