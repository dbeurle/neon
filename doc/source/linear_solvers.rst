Linear Solvers
==============

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

To specify an iterative solver require additional fields due to the white-box nature of the methods.  If these are not set, then defaults will be chosen for you.  The following table demonstrates the defaults

.. table:: Iterative solvers defaults
   :widths: auto

   ======================== ============================================
   Iterative solver option  Default value
   ======================== ============================================
   ``"preconditioner"``     ``"diagonal"``
   ``"backend"``            ``"cpu"``
   ``"tolerance"``          ``1.0e-6``
   ``"maximum_iterations"`` ``1.0e-6``
   ======================== ============================================

and the option meanings

.. table:: Iterative solvers available ``"Type" : "keyword"``
   :widths: auto

   ==================== ============================================
   Additional options   Details
   ==================== ============================================
   ``"backend"``        ``"cpu"`` and ``"gpu"`` for selecting computation backend (``"gpu"`` requires `-DENABLE_OCL=1` or `-DENABLE_CUDA` during compile time)
   ``"preconditioner"`` ``"diagonal"`` is the only supported preconditioner available
   ==================== ============================================

An example of an iterative solver definition ::

     "linear_solver" {
         "type" : "iterative",
         "preconditioner" : "diagonal",
         "backend" : "gpu",
         "tolerance" : 1.0e-6,
         "maximum_iterations" : 1500
     }

All linear solvers use double floating point precision which may incur performance penalties on GPU devices.
