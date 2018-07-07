[![Build Status](https://travis-ci.org/dbeurle/neon.svg?branch=master)](https://travis-ci.org/dbeurle/neon)
[![Coverage Status](https://coveralls.io/repos/github/dbeurle/neon/badge.svg?branch=master)](https://coveralls.io/github/dbeurle/neon?branch=master)
[![Documentation Status](https://readthedocs.org/projects/neon-user/badge/?version=latest)](http://neon-user.readthedocs.io/en/latest/?badge=latest)
[![Documentation](https://codedocs.xyz/dbeurle/neon.svg)](https://codedocs.xyz/dbeurle/neon/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# neon
A non-linear finite element code.  This project is still under development and is not considered stable.

## Building

### Docker

Go to the top level directory in the git repository and enter

- `sudo docker build -t neon .`

This will pull down the base image and then clone the git repository inside the docker container.  The original hosting can be found at https://hub.docker.com/r/dbeurle/neon/

Once this is completed, enter the following commands to enter the docker container

- `sudo docker run -i -t neon /bin/bash`

The instructions in the section below can now be followed with the copy of the repository inside the docker container without worrying about dependencies.

### Linux

The external dependencies are:
 - Boost filesystem (to be replaced with `std::filesystem`)
 - Pastix and MUMPS for direct linear solvers
 - VTK for writing out simulation data to a ParaView compatible format
 - An OpenMP enabled C++17 compiler
 - Intel Thread Building Blocks TBB

Other dependencies are pulled in during build time with `CMake` and include

 - Eigen
 - range-v3
 - Termcolor for colour terminal support
 - Catch for unit testing

For best performance build with native optimisations on for the machine you are using.  This will automatically trigger Eigen to use wider SIMD instruction when available on native architecture at the expense of portability.

The build instructions for debug are
- `$ mkdir build && cd build/`
- `$ cmake ..`
- `$ make all`

and for release mode with optimisations

- `$ mkdir build && cd build/`
- `$ cmake -DCMAKE_BUILD_TYPE=Release ..`
- `$ make all -j<# cpus>`
- `$ sudo make install`

If the code is only used on a particular machine, then machine specification optimisations can be invoked by specifying the `CMake` symbol
- `-DENABLE_NATIVE=1`

which enables the compiler to generate the best machine code for your platform.

If you have an NVIDIA graphics card, then you can use the CUDA enabled iterative solvers by specifying the `CMake` symbol

- `-DENABLE_CUDA=1`

For checking the successful compilation of the program, invoke the test suite by executing

- `ctest`

in the build directory.

#### Ubuntu 18.04

Install dependencies through the package manager:

`$ sudo apt install cmake git mercurical libz-dev libcurl3-dev libvtk6-dev libtbb-dev libboost-filesystem-dev libmumps-seq-dev libopenblas-dev libarpack2-dev, libscotch-dev, hwloc libhwloc-dev libgfortran-7-dev`

Then clone the repository and add

`$ git clone <neon>`

and enter the repository

`$ cd neon/docker-base && sh install_pastix.sh`

After this compiles and fails to install, enter the commands to install and link the libraries

`cd pastix_5.2.3/build && sudo make install && sudo ln -s /usr/local/lib/libpastix.so /usr/lib/libpastix.so`

Provide the `sudo` password when prompted.  Go back to the top `neon` directory and create and enter the `build` directory

`$ mkdir build && cd build`

let `CMake` know want to use `g++-7` as a C++ compiler

`$ export CXX=g++-7; cmake ..`

finally compiling with

`$ make all -jN`

where `N` is the number of parallel build jobs you want to run.

## Licensing

See the LICENSE.md file for the project license and the licenses of the included dependencies.
