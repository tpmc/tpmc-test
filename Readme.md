# Tests for the topology preserving marching cubes (TPMC)

[![Build Status](https://travis-ci.org/tpmc/tpmc-test.svg?branch=master)](https://travis-ci.org/tpmc/tpmc-test)

The test programs are build using Cmake.

In addition to the TPMC library, the test programs need the following third party libraries:

- Eigen3
- Boost

Create a build directory and call cmake:

```
mkdir build;
cd build;
cmake ..
```

When the cmake call was successfull, you can build and run the test problems by calling

```
make check;
```

in the build directory. The result of the test is stored in a file `test<n>.log` where `<n>` denotes the test number. If you want to compile the test programs without running them, call

```
make;
```

The tests are located in the `src` subdirectory and can be run by calling

```
./<test-name> output.log;
```

where `<test-name>` denotes the name of the test executable.
