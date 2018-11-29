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

Eigenvalue problems
-------------------

The solution of an eigenvalue problem arises in natural frequency analysis and buckling analysis in structural dynamics.  Usually only the solution of a few eigenvalues are required and this becomes computationally feasible for large sparse matrices.  These eigenvalues and eigenvectors correspond to the natural frequency and the mode shape of the structure for natural frequency analysis.  For the linear elastic buckling load of a structure, the eigenvalues provide the limit load and the deformed shape.

Unlike direct methods for linear systems, eigenvalues have to be solved for in an iterative fashion.  There are several different algorithms available depending on the problem.

An example of an eigenvalue solver definition ::

     "eigen_solver" {
         "type" : "lanczos",
         "values" : 15,
         "spectrum" : "lower"
     }

where the ``"type"`` field indicates what algorithm to use (``"arpack"``, ``"lanczos"`` and ``"power_iteration"``) are available.  The ``"values"`` keyword determines how many eigenvalues are to be solved for.  This should be much less than the total number of degrees of freedom in the system.  Finally the ``"spectrum"`` keyword indicates from which end of the spectrum the values will be computed, where ``"lower"`` indicates Eigenvalues from the lowest frequency and ``"upper"`` computes the higher frequency Eigenvalues.
