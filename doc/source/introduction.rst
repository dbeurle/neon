Introduction
============

neon is a finite element solver written in C++.  It aims to solve industrial style finite element problems exploiting modern multicore architectures such as GPUs to accelerate simulations.  The finite element problem definition is specified through a ``json`` format input file.  Information regarding material, parts, constitutive models, boundary conditions and solution techniques are specified in this input file.

The mesh file is also given in ``json`` format which makes it trivial to parse and provides a standardised mesh format that includes syntax highlighting and code folding features.

A containerised version of neon is also available at Dockerhub, which allows the build environment to be simply reproduced.
