# Tests for the topology preserving marching cubes (TPMC)

[![Build Status](https://travis-ci.org/tpmc/tpmc-test.svg?branch=master)](https://travis-ci.org/tpmc/tpmc-test)

The test programs are build using Cmake.

The test programs need the following third party libraries:

- Eigen3
- Boost
- python-numpy

To get started, clone the git repository, create a build directory and call cmake:

```
git clone https://github.com/tpmc/tpmc-test.git <your-tpmc-test-dir>
cd <your-tpmc-test-dir>
mkdir build
cd build
cmake ..
```

When the cmake call was successfull, you can build and run the test problems by calling

```
make check
```

in the build directory. The result of the test is stored in a file `test<n>.log` where `<n>` denotes the test number. If you want to compile the test programs without running them, call

```
make
```

The tests are located in the `src` subdirectory and can be run by calling

```
./<test-name> output.log
```

where `<test-name>` denotes the name of the test executable.

License
-------

The TPMC library, headers and test programs are free open-source
software, dual-licensed under version 3 or later of the GNU Lesser
General Public License and version 2 of the GNU General Public License
with a special run-time exception.

See the file LICENSE.txt for full copying permissions.
