FROM fedora:28

RUN dnf update -y
# Utilities
RUN dnf install -y wget bzip2 patch git mercurial cmake lcov
# Compilers
RUN dnf install -y clang clang-tools-extra libomp-devel gcc-c++
# Scientific libraries
RUN dnf install -y openblas-devel MUMPS-devel boost-devel vtk-devel hwloc-devel scotch-devel tbb-devel arpack-devel catch-devel viennacl-devel
# Compile and install PaStiX 5
COPY docker/add_libraries.diff docker/cmake_patch.diff docker/install_pastix.sh /
RUN sh install_pastix.sh
# OpenCL vendor loader
RUN dnf install -y ocl-icd-devel
# OpenCL runtimes for xeon-based build server
RUN wget http://registrationcenter-download.intel.com/akdlm/irc_nas/12556/opencl_runtime_16.1.2_x64_rh_6.4.0.37.tgz

RUN tar -xf opencl_runtime_16.1.2_x64_rh_6.4.0.37.tgz
COPY docker/ocl_config.diff opencl_runtime_16.1.2_x64_rh_6.4.0.37
RUN cd opencl_runtime_16.1.2_x64_rh_6.4.0.37/ && patch -p0 < ocl_config.diff && sh install.sh --silent silent.cfg

# Coverage scripts
RUN pip3 install gcovr cpp-coveralls && pip install pyyaml
