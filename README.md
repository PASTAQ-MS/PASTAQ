# Test PASTAQ on `binder`

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PASTAQ-MS/PASTAQ/HEAD?labpath=examples%2FpastaqExamples.ipynb)

# Installation

## Python 3 virtual environment

```sh
pip install pastaq
```

## Installing from source

You need to install a suitable C++ compiler and
corresponding build tools for your platform as well as CMake. 
The instructions listed here refer to the installation of PASTAQ's Python
bindings.  Currently the only external dependencies, including zlib, are included as git submodules.  

To get started, clone this repository and initialize git submodules:

```
git clone https://github.com/PASTAQ-MS/PASTAQ.git
cd PASTAQ
git submodule init
git submodule update --remote
```

As usual, it is strongly recommended to create a **Python 3 environment** in which to build Pastaq, and the core development has been with Python `3.9`, but `3.10`, `3.11` and `3.12` should also work.

```sh
python -m pip install --upgrade pip
python -m pip install build
python -m pip install wheel

# create the .whl file in the ./dist folder
python -m build --installer pip --wheel
```

## Windows

When building Pastaq in Windows, it may be helpful to first open a Visual Studio command prompt using Tools->Visual Studio Command Prompt in the Visual Studio IDE so that you have access to the compiler and linker.  Then, in that command window, activate your PASTAQ Python environment and proceed with the instructions.

#### Powershell

```
Get-ChildItem ./dist/*.whl | ForEach-Object { pip install $_.FullName }
```

#### CMD command prompt

```
for %f in (./dist\*.whl) do pip install %f
```
## Linux

```
find ./dist/*.whl | xargs pip install 
```

<!--- # # Build the module and install it in your system: -->
<!--- ```sh -->
<!--- # Installation -->
<!--- python3 setup.py install --user -->

<!---  # Development -->
<!--- python3 setup.py develop --user -->
<!--- ``` -->

Now it can be imported and used in python as follows:

```python
import pastaq
raw_data = pastaq.read_mzxml(...)
```

# Usage

Examples of the usage of the PASTAQ can be found in the `examples` folder. To
run them, install pastaq as previously described, update the input path of the
mzXML and mzID files, change any necessary parameters and run it with:

```
python examples/small_range.py
```

You can use any mzXML files and identifications in mzIdentML v1.1+. If no
identifications are available, remove the `ident_path` from the input files
array or set it to `'none'`. You can find the files we used for testing and
development via ProteomeXchange, with identifier [PXD024584](https://www.ebi.ac.uk/pride/archive/projects/PXD024584).

Processing of mzML files is in an early stage and may lead to some issues.

For more information about PASTAQ and the configuration of the parameters,
please visit [the official website][website].

[website]: https://pastaq.horvatovichlab.com

<!-- # Compile the pastaq library from source -->

<!-- For an out of source build of the library, create a build directory and run cmake: -->

<!-- ```sh -->
<!-- mkdir build -->
<!-- cd build -->
<!-- cmake .. -->
<!-- make -->
<!-- ``` -->

<!-- If you wish to enable the compilation of the tests you need to set up the -->
<!-- `PASTAQ_ENABLE_TESTS` flag to 1. A limited number of automated test are -->
<!-- currently available but we are looking forward to improve this in the future. -->

<!-- ```sh -->
<!-- mkdir build -->
<!-- cd build -->
<!-- cmake .. -DPASTAQ_ENABLE_TESTS=1 -->
<!-- make -->
<!-- make test -->
<!-- ``` -->

<!-- Additionally, you can use the Ninja building tool for faster compilation times. -->

<!-- ```sh -->
<!-- mkdir build -->
<!-- cd build -->
<!-- cmake .. -DPASTAQ_ENABLE_TESTS=1 -GNinja -->
<!-- ninja -->
<!-- ninja test -->
<!-- ``` -->

## Third-Party Libraries

This project uses the TDF SDK (v2.21.0.x) to interface with timsTOF data. The `include` folder contains essential headers required for compilation.

### File Structure
- `libtimsdata.so` (Linux) or `timsdata.dll` (Windows): Shared libraries.
- `include/`: Header files required to build and use the SDK.
  - `c/`: Contains C/C++ header files.
    - `timsdata.h`: Main header for interfacing with timsTOF data.
    - `timsvisualization.h`: Header for visualization functionalities.
    - Other headers as required by the SDK.

These files are distributed in compliance with the SDK's licensing terms.

# Tutorial Jupyter notebook

Tutorial demonstrating PASTAQ functionality to read and process LC-MS/MS data is available in examples/pastaqExamples.ipynb. Click the Binder icon above to try out the tutorial in Binder virtual environment.

# How to cite this work

The main manuscript has been published in as Open Access Analytical Chemistry with the following details: [Alejandro Sánchez Brotons, Jonatan O. Eriksson, Marcel Kwiatkowski, Justina C. Wolters, Ido P. Kema, Andrei Barcaru, Folkert Kuipers, Stephan J. L. Bakker, Rainer Bischoff, Frank Suits, and Péter Horvatovich, Pipelines and Systems for Threshold-Avoiding Quantification of LC–MS/MS Data, Analytical Chemistry, 2021, 93, 32, 11215–11224](https://pubs.acs.org/doi/10.1021/acs.analchem.1c01892).
