
FROM fedora:28

# Install any needed packages specified in dependencies.txt
RUN dnf update -y

RUN dnf install -y  wget bzip2 patch git mercurial cmake \
                    gcc-c++ openblas-devel MUMPS-devel \
                    boost-devel vtk-devel hwloc-devel scotch-devel \
                    tbb-devel arpack-devel gcovr

# Copy the current directory contents into the container at /
COPY docker/add_libraries.diff docker/cmake_patch.diff docker/install_pastix.sh /

RUN sh install_pastix.sh
