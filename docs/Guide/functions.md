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

---

