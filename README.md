# neon
A non-linear finite element code.  This project is still under development and is not considered stable.

## Building

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
- `make all`

For checking the successful compilation of the program, invoke the test suite by

`ctest`

## Licensing

See the LICENSE.md file for the project license and the licenses of the included dependencies.
