# Installation

Clone this repository and initialize git submodules:

```
git clone git@git.horvatovichlab.com:PASTAQ/pastaq.git
git submodule init
git submodule update
```

Build the module and install it in your system:

```sh
# Installation
python3 setup.py install --user

# Development
python3 setup.py develop --user
```

Now it can be imported and used in python as follows:

```python
import pastaq
raw_data = pastaq.read_mzxml(...)
```

# Compile the pastaq library from source

For an out of source build of the library, create a build directory and run cmake:

```sh
mkdir build
cd build
cmake ..
make
```

If you wish to enable the compilation of the tests you need to set up the
`PASTAQ_ENABLE_TESTS` flag to 1.

```sh
mkdir build
cd build
cmake .. -DPASTAQ_ENABLE_TESTS=1
make
make test
```

You can use the Ninja building tool for faster compilation times.

```sh
mkdir build
cd build
cmake .. -DPASTAQ_ENABLE_TESTS=1 -GNinja
ninja
ninja test
```
