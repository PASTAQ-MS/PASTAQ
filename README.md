# Installation

Build the module for development with:

```sh
python3 setup.py develop --user
```

Now it can be imported and used in python as follows:

```python
import tapp
raw_data = tapp.read_mzxml(...)
```

# Compile from source

For an out of source build of the library, create a build directory and run cmake:

```sh
mkdir build
cd build
cmake ..
make
```

If you wish to enable the compilation of the tests you need to set up the
`TAPP_ENABLE_TESTS` flag to 1.

```sh
mkdir build
cd build
cmake .. -DTAPP_ENABLE_TESTS=1
make
make test
```

You can use the Ninja building tool for faster compilation times.

```sh
mkdir build
cd build
cmake .. -DTAPP_ENABLE_TESTS=1 -GNinja
ninja
ninja test
```
