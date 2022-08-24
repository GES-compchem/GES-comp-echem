# User guide

The `GES-comp-echem` library is constructed based on the following submodules:

* `systems`: classes for storing informations about the system under study, such as geometry, number of atoms, calculated energies, charge, spin, etc.

* `tools`: various useful functions called internally throughout the code. The user should generally not need to use this module during normal operation.

* `wrappers`: wrappers for interfacing with external software, such as **Orca**, **DFTB+**, **CREST**, etc.

* `functions`: mid-level algorithms for manually calculating specific physical observables, such as pKa or redox potentials. The user must provide all the necessary species for the calculation, for example the specific protonated and deprotonated species in a pKa calculation.

* `algorithms`: high-level algorithms for automatically calculating specific physical observables, considering all the "non-trivial" complications of the case. For example, taking into account protonation states and PCET mechanisms in the calculation of the redox potential. The user only needs to provide a "generic" state of the system of interest, and the algorithm should find the exact species involved in the calculation (e.g., the correct tautomer/conformer at the correct protonation state at a given pH)
