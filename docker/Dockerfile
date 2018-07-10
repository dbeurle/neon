
FROM fedora:28

# Install any needed packages specified in dependencies.txt
RUN dnf update -y
# Utilities
RUN dnf install -y wget bzip2 patch git mercurial cmake lcov
# Compilers
RUN dnf install -y clang clang-tools-extra libomp-devel gcc-c++
# Scientific libraries
RUN dnf install -y openblas-devel MUMPS-devel boost-devel vtk-devel hwloc-devel scotch-devel tbb-devel arpack-devel
# Coverage scripts
RUN pip3 install gcovr cpp-coveralls

# Copy the current directory contents into the container at root directory /
COPY docker/add_libraries.diff docker/cmake_patch.diff docker/install_pastix.sh /
RUN sh install_pastix.sh
