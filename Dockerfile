# Use the base image to avoid multiple dependency installations
FROM dbeurle/neon:base

# Pull down the code
RUN git clone --depth=1 https://github.com/dbeurle/neon.git
