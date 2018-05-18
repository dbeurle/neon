
FROM fedora:latest

# Install any needed packages specified in dependencies.txt
RUN dnf update -y

RUN dnf install -y  wget bzip2 patch lcov git mercurial cmake \
                    gcc-c++ openblas-devel MUMPS-devel \
                    boost-devel vtk-devel hwloc-devel scotch-devel \
                    tbb-devel arpack-devel

# Set the working directory to /
WORKDIR /

# Copy the current directory contents into the container at /
ADD . /

RUN sh install_pastix.sh
