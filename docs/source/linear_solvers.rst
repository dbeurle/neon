Linear solvers
==============

Linear solvers make up the largest part of the solution time in the finite element and we provide many options to the user.  Large scale systems should use iterative solvers as their memory usage is very small in comparison to their direct solver brethren.

Neon uses a selection of linear solvers, most of which are multithreaded.  Direct solvers require no options apart from the solver name.

.. table:: Direct solvers available ``"Type" : "keyword"``
   :widths: auto

   ============ ============================================
   Type keyword Details
   ============ ============================================
   Direct       Inbuilt from the Eigen library (LLT or LU)
   PaStiX       Interfaces to the PaStiX library (LLT or LU)
   MUMPS        Interfaces to the MUMPS library (LDLT or LU)
   ============ ============================================

On the other hand iterative solvers require the maximum number of iterations and the maximum tolerance.  If these are not set, then defaults will be chosen for you.  Iterative methods specify the following fields ::

    "LinearSolver" {
        "Type" : "Iterative",
        "Tolerance" : 1.0e-6,
        "MaxIterations" : 1500
    }

.. table:: Iterative solvers available ``"Type" : "keyword"``
   :widths: auto

   ============ ============================================
   Type keyword Details
   ============ ============================================
   Iterative    Inbuilt from the Eigen library (Diagonal preconditioned CG or BiCGStab)
   IterativeGPU GPU accelerated diagonal preconditioned CG (requires -DENABLE_CUDA=1 in build)
   ============ ============================================
