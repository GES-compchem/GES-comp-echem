# Contributor's guide

`GES-comp-echem` is an open source project and, as such, all the contributions are well accepted. If you want to contribute to the library, a simple guide on how to interface with the GitHub repository is provided in what follows.

## General development process

* If you are a first-time contributor:

    * Go to [https://github.com/GES-compchem/GES-comp-echem](https://github.com/GES-compchem/GES-comp-echem) and click the “fork” button to create your own copy of the project repository.

    * Clone the project to your local computer using the command (where `<YOUR_USERNAME>` represent your personal GitHub username): 
        ```
        git clone https://github.com/<YOUR_USERNAME>/GES-comp-echem
        ```

    * Enter the reposiotry directory using the command:
        ```
        cd GES-comp-echem
        ```

    * Add the upstream repository using the command:
        ```
        git remote add upstream https://github.com/GES-compchem/GES-comp-echem.git
        ```

    * Now, when running the `git remote -v` command, the following reposiotries shold be visible:
        * `upstream`, which refers to the GES-compchem repository
        * `origin`, which refers to your personal fork

* Developing your contributions:
    
    * Pull the latest changes from upstream:
        ```
        git checkout main
        git pull upstream main
        ```

    * Create a branch for the feature you want to work on. Since the branch name will appear in the merge message, use a sensible name.:
        ```
        git checkout -b the_name_of_the_branch
        ```

    * Commit locally as you progress (`git add` and `git commit`) Use a properly formatted commit message. If possible, write tests that fail before your change and pass afterward, run all the tests locally. Be aware that the whole suite of tests can be run using `pytest --cov`. Before the commit use `tox` to verify the compatibility with all the supported version of python (`tox` will run only unit tests if executed using the provided `tox.ini` configuration file). Be sure to document any changed behavior in docstrings, keeping to the [NumPy docstring standard](https://numpydoc.readthedocs.io/en/latest/format.html). Sphinx annotation are very welcome. If brand new functionality are added to the package make sure to update the documentation as well.

* To submit your contribution:

    * Push your changes back to your fork on GitHub:
        ```
        git push origin the_name_of_the_branch
        ```
    * Follow the GitHub authentication process.

    * Go to GitHub. The new branch will show up with a green Pull Request button. Make sure the title and message are clear, concise, and self- explanatory. Then click the button to submit it.

    * Ask for review from the development team.

# Codebase overview

The purpose of the library is to give to the user a simple, self-contained, interface to computational chemistry software, tools and internal routines useful in carrying out computational chemistry calculations. 

The library is built around the `System` object. This object, defined in `compechem.systems`, represents a molecule that is defined by a structure and a set of properties.

The library articulate all its functions across different modules each of which carrying different functionality:

* The `core` module contains all the basic components required to operate the library and define its objects. A detailed description of the inner working of the `core` module can be found in the [corresponding API page](API-core).

* The `engine` module contains all the interface to computational chemistry softares that, once given a molecule, are capable of computing and setting `System` properties. A detailed description of the inner working of the `engine` module can be found in the [corresponding API page](API-engines).

* The `wrapper` module contains all the interface to all the computational chemistry softares that does not match the `engine` module definition. As an example, wrappers to softwares like `crest` and `packmol` that operates or generate `System` object without computing properties can be found here. A detailed description of the inner working of the `wrapper` module can be found in the [corresponding API page](API-wrappers).

* The `functions` module contains all the functions encoding workflows, typically originated combining operation performed by the engines, to "manually" calculate some physical properties of the system under study. As an example the computation of the pKa of a system or its Fukui functions can be found here. A detailed description of the inner working of the `functions` module can be found in the [corresponding API page](API-functions).

* The `tools` module contains all the auxiliary tools used to perform specific operations like parsing files, converting file formats, providing an interface to visualization tools etc. A detailed description of the inner working of the `tools` module and all its components can be found in the [corresponding API page](API-tools).

# Directory Structure

The directroy structure of the repository is the following:

```
.
├── bld.dat                 \\ File required by the conda package manager
├── build.sh                \\ File required by the conda package manager
├── compechem               \\ Folder containing the package code with all the modules
├── docs                    \\ Folder containing the files required to build the jupyter-books documentation 
├── LICENSE                 \\ LICENSE file
├── meta.yaml               \\ Configuration file defining the conda package properties and dependencies
├── pyproject.toml          \\ Configuration file for the build system
├── README.md               \\ README file
├── requirements_dev.txt    \\ Requirements for testing and development
├── requirements.txt        \\ Requirements for building and using the package
├── setup.cfg               \\ Setup configuration file
├── setup.py                \\ The setup script
├── tests                   \\ The folder containg the packages tests
├── tox.ini                 \\ The tox configuration file to be used for automatic testing
└── validation              \\ Legacy folder used in the past to build example and test references
```

The package folder is organized as follows:

```
.
├── config.py                   // Configuration file setting the library configuration (MPI flags, opeation mode, ecc..)
├── constants.py                // Script containing all the physical constants and element data
├── core
│   ├── base.py                 // The definition of the Engine base class used to derive the engine modules
│   ├── dependency_finder.py    // The script used to locate all the third party softwares required
│   ├── geometry.py             // The module encoding the geometric structure of a given molecule
│   ├── __init__.py
│   ├── properties.py           // The module containing and organizing the properties associated to a given System
│   └── spectroscopy.py         // The module containing all the information related to the spectroscopic response of a system
├── engines
│   ├── dftbplus.py             // The module containing the definition of the DFTBInput engine
│   ├── __init__.py 
│   ├── orca.py                 // The module containing the definition of the OrcaInput engine
│   └── xtb.py                  // The module containing the definition of the XtbInput engine
├── functions
│   ├── fukui.py                // The module containing all the functions used to compute Fukui functions
│   ├── __init__.py
│   ├── pka.py                  // The module containing all the functions capable of computing the pKa of a given molecule
│   └── potential.py            // The module containing the functions capable of computing the reduction potential of a given system 
├── __init__.py
├── systems.py                  // The module where the System and Ensamble objects are defined
├── tools
│   ├── cubetools.py            // The module containing the definition of a Cube object
│   ├── externalutilities.py    // The module containing a miscellaneous of different tools to manipulate files
│   ├── __init__.py
│   ├── internaltools.py        // The module containing a miscellaneous of different internal tools
│   ├── moglitools.py           // The interface to the Mogli molecular viewer package
│   ├── reorderenergies.py      // Tool used to reorder systems based on their energy
│   ├── vmdtools.py             // The interface to the VMD renderer
│   └── xyz2mol.py              // 
└── wrappers
    ├── crest.py                // The module containing the definition of the CREST wrapper
    ├── __init__.py
    └── packmol.py              // The module containing the definition of the PACKMOL wrapper
```

The detailed information about the working of each module is provided in the [API documentation](API-main).