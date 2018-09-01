FROM fedora:28

RUN dnf update -y
# Utilities
RUN dnf install -y wget bzip2 patch git mercurial cmake lcov
# Compilers
RUN dnf install -y clang clang-tools-extra libomp-devel gcc-c++
# OpenCL runtimes
RUN dnf install -y pocl-devel
# Scientific libraries
RUN dnf install -y openblas-devel MUMPS-devel boost-devel vtk-devel hwloc-devel scotch-devel tbb-devel arpack-devel catch-devel viennacl-devel
# Coverage scripts
RUN pip3 install gcovr cpp-coveralls
RUN pip install pyyaml

# Compile and install PaStiX 5
COPY docker/add_libraries.diff docker/cmake_patch.diff docker/install_pastix.sh /
RUN sh install_pastix.sh
