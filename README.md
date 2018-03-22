[![Build Status](https://travis-ci.org/dbeurle/neon.svg?branch=master)](https://travis-ci.org/dbeurle/neon)
[![Coverage Status](https://coveralls.io/repos/github/dbeurle/neon/badge.svg?branch=master)](https://coveralls.io/github/dbeurle/neon?branch=master)
[![Documentation Status](https://readthedocs.org/projects/neon-user/badge/?version=latest)](http://neon-user.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# neon
A non-linear finite element code.  This project is still under development and is not considered stable.

## Building

### Docker build

Go to the top level directory in the git repository and enter

- `sudo docker build -t neon .`

This will pull down the base image and then clone the git repository inside the docker container.  The original hosting can be found at https://hub.docker.com/r/dbeurle/neon/

Once this is completed, enter the following commands to enter the docker container

- `sudo docker run -i -t neon /bin/bash`

The instructions in the section below can now be followed with the copy of the repository inside the docker container without worrying about dependencies.

### Linux

The external dependencies are:
 - Boost filesystem
 - Pastix and MUMPS for direct linear solvers
 - Scotch for matrix reordering
 - OpenBLAS for linking to linear solvers
 - VTK for writing out simulation data to ParaView
 - A c++17 compatible compiler with unicode support (clang 4.0 or clang 5.0 development)
 - Jsoncpp for processing mesh files
 - OpenMP for parallelisation

Other dependencies are pulled in during build time with `CMake` and include

 - Eigen
 - Range v3
 - Termcolor for colour terminal support
 - Catch for unit testing

For best performance build with native optimisations on for the machine you are using.  This will automatically trigger Eigen to use wider SIMD instruction when available on native architecture at the expense of portability.

The build instructions for debug are
- `mkdir build && cd build/`
- `cmake ..`
- `make all`

and for release mode with optimisations

- `mkdir build && cd build/`
- `cmake -DCMAKE_BUILD_TYPE=Release ..`
- `make all -j<# cpus>`
- `make install`

If the code is only used on a particular machine, then machine specification optimisations can be invoked by specifying the `CMake` symbol
- `-DENABLE_NATIVE=1`

which enables the compiler to generate the best machine code for your platform.

If you have an NVIDIA graphics card, then you can use the CUDA enabled iterative solvers by specifying the `CMake` symbol

- `-DENABLE_CUDA=1`

For checking the successful compilation of the program, invoke the test suite by executing

- `ctest`

in the build directory.



## Licensing

See the LICENSE.md file for the project license and the licenses of the included dependencies.
