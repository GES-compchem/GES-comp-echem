The `GES-comp-echem` package can be installed in a Conda environment with the command:
```
conda install -c greenenergystorage GES-comp-echem
```

The library can be imported in a Python script as a whole via the following syntax:

```python
import compechem
```

Alternatively, individual submodules/classes/functions can be imported separately:

```python
from compechem import systems
from compechem.wrappers import dftbplus
from compechem.wrappers.packmol import packmol_cube
```

For a more detailed explanation of the available features in each submodule, please refer to their specific page in this User Guide.