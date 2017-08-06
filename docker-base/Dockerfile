# Use Fedora 26 as an official parent image
FROM fedora:26

# Install any needed packages specified in dependencies.txt
RUN dnf install -y  wget bzip2 patch lcov git mercurial cmake \
                    gcc-c++ openblas-devel MUMPS-devel jsoncpp-devel \
                    boost-devel vtk-devel hwloc-devel scotch-devel

# Set the working directory to /
WORKDIR /

# Copy the current directory contents into the container at /
ADD . /

RUN sh install_pastix.sh
