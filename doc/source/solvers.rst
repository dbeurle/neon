Solvers
=======

Linear systems
--------------

Linear solvers typically take up the largest portion of the solution time in the finite element method.  Very large scale systems of million degrees of freedom should use iterative solvers as their memory usage is small in comparison to their direct solver brethren.  For problems where it is possible to use a direct solver then this is the suggested default.

Neon uses a selection of linear solvers, most of which employ shared memory parallelisation or have GPU-accelerated counterparts.  Direct solvers require no options apart from the solver name.  The type of solver (symmetric or unsymmetric) is automatically deduced based on boundary conditions, constitutive model and other algorithmic choices.

.. table:: Direct solvers available ``"type" : "keyword"``
   :widths: auto

   ============ ============================================
   Type keyword Details
   ============ ============================================
   ``"direct"`` Inbuilt from the Eigen library (LLT or LU)
   ``"PaStiX"`` Interfaces to the PaStiX library (LLT or LU)
   ``"MUMPS"``  Interfaces to the MUMPS library (LDLT or LU)
   ============ ============================================

An example of using the ``"PaStiX"`` solver is ::

    "linear_solver" {
        "type" : "PaStiX"
    }

To specify an iterative solver require additional fields due to the white-box nature of the methods.  If these are not set, then defaults will be chosen for you.  The following table demonstrates the defaults, where each iterative solver currently uses a diagonal pre-conditioner due to ease of parallelisation

.. table:: Iterative solvers defaults
   :widths: auto

   ======================== ============================================
   Iterative solver option  Default value
   ======================== ============================================
   ``"device"``             ``"cpu"``
   ``"tolerance"``          ``1.0e-6`` (Relative tolerance)
   ``"maximum_iterations"`` ``1.0e-6`` (Maximum iterations before failure)
   ======================== ============================================

and the option meanings

.. table:: Iterative solvers available ``"Type" : "keyword"``
   :widths: auto

   ==================== ============================================
   Additional options   Details
   ==================== ============================================
   ``"device"``         ``"cpu"`` and ``"gpu"``
   ``"backend"``        Provide options for the acceleration framework
   ==================== ============================================

Note than when selecting ``"gpu"`` the binary must have compiled in support that is enabled through the ``-DENABLE_OPENCL=1`` or ``-DENABLE_CUDA`` CMake flags.

An example of an iterative solver definition using the CUDA iterative linear solver ::

     "linear_solver" {
         "type" : "iterative",
         "device" : "gpu",
         "tolerance" : 1.0e-6,
         "maximum_iterations" : 1500,
         "backend" : {
            "type" : "cuda",
            "device" : 0
         }
     }

The ``"type"`` field can also contain ``"opencl"`` to use the OpenCL kernels for computation.  Because of the multi-device nature of OpenCL a ``"device"`` and a ``"platform"`` are then required to select the desired device ::

    "linear_solver" {
        "type" : "iterative",
        "device" : "gpu",
        "tolerance" : 1.0e-6,
        "maximum_iterations" : 1500,
        "backend" : {
           "type" : "opencl",
           "platform" : 0,
           "device" : 0
        }
    }

All linear solvers use double floating point precision which may incur performance penalties on GPU devices however testing shows a significant increase in performance is still obtained with double precision due to the higher memory bandwidth.  A single precision version of Krylov subspace solvers is not recommended due to round-off error in computing the search direction.

Eigenvalue problems
-------------------

The solution of an eigenvalue problem arises in natural frequency analysis in structural dynamics.  Usually only the solution of a few eigenvalues are required and is computationally feasible.  These eigenvalues and eigenvectors correspond to the natural frequency and the mode shape of the structure for natural frequency analysis.  For the linear elastic buckling load of a structure, the eigenvalues provide the limit load and the deformed shape.

Unlike direct methods for linear systems, eigenvalues have to be solved for in an iterative fashion.  There are several different algorithms available depending on the problem.

An example of an eigenvalue solver definition ::

     "eigenvalue_solver" {
         "type" : "lanczos",
         "values" : 15,
         "sort" : "lowest"
     }

where the ``"type"`` field indicates what algorithm to use and ``"values"`` determines how many eigenvalues are to be solved for.  This should be much less than the total number of degrees of freedom in the system.
