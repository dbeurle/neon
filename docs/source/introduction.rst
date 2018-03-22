Introduction
============

neon is a finite element solver written in C++.  It aims to solve industrial style finite element problems exploiting modern multicore architectures such as GPUs for accelerating simulations.  The finite element problem definition is specified through a ``json`` format input file.  Here material, parts, constitutive models, boundary conditions and solver information is specified.

The mesh file is also given in ``json`` format which makes it trivial to parse and provides a standardised mesh format that includes syntax highlighting and code folding features.

A containerised version of neon is also available at Dockerhub, so you can easily reproduce the build environment.
