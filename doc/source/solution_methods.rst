****************
Solution Methods
****************

Linear Equilibrium and Transient
================================

Solving linear problems involves one invocation of a linear solver and one assembly step and is therefore inexpensive to perform.  These routines are automatically selected based on the problem and deserve no special discussion.

Non-linear Equilibrium
======================

In the non-linear finite element framework, the solution of the non-linear equations requires linearisation and an incremental approach to determine the solution.  This means that multiple invocations of a linear solver are required to reach a specified equilibrium.

To perform the solution stage, neon implements the full Newton-Raphson method such that an updated tangent matrix is computed in each iteration.  Since the assembly of the tangent stiffness matrix is implemented in parallel, the computational cost is very low in comparison to the cost of a linear solve.  Other finite element solvers will avoid the computation of stiffness matrix due to the computational cost at the expense of improved convergence properties.

The iterative nature of a non-linear problem requires the use of tolerances to determine if the results are sufficiently converged.  For this, non-linear simulation cases need to specify the relative displacement, force residuals and the maximum number of Newton-Raphson iterations to perform before a cutback.  For additional control, the relative or absolute tolerances can be chosen based on the physics of the problem ::

    "nonlinear_options" : {
        "displacement_tolerance" : 1.0e-5,
        "residual_tolerance" : 1.0e-4,
        "absolute_tolerance" : true,
        "newton_raphson_iterations" : 15
    }

Methods to improve the properties of the Newton-Raphson could be implemented on top of the current non-linear solvers, such as line searching algorithms to improve convergence properties.


Non-linear Implicit Dynamic
===========================

The code to solve the implicit dynamic problem has been implemented but not yet tested.  If required, please provide a test case to the developers.

Natural Frequency
=================

The generalised Eigenvalue problem can be solved by specifying ::

    "steps" : [{
        ...
        "solution" : "natural_frequency",
        ...
    }]

Note that loads do not need to be specified for free vibration analysis but the restraints must be specified.  Since this solution method builds a mass matrix, the material for the mesh must specify a material density.

For a selection of valid Eigenvalue solvers, please refer to the solvers section of this manual.
