# Use the base image to avoid multiple dependency installations
FROM neon:base

# Pull down the code
RUN git clone --depth=1 https://github.com/dbeurle/neon.git
RUN cd neon && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make all -j4
