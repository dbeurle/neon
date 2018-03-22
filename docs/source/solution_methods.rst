****************
Solution methods
****************

Linear equilibrium and transient
================================

Solving linear problems involves one invocation of a linear solver and one assembly step and is therefore inexpensive to perform.  These routines are automatically selected based on the problem and deserve no special discussion.

Non-linear equilibrium
======================

In the non-linear finite element framework, the solution of the non-linear equations requires linearisation and an incremental approach to determine the solution.  This means that multiple invocations of a linear solver are required to reach a specified equilibrium.

To perform the solution stage, neon implements the full Newton-Raphson method such that an updated tangent matrix is computed in each iteration.  Since the assembly of the tangent stiffness matrix is implemented in parallel, the computational cost is very low in comparison to the cost of a linear solve.  Other finite element solvers will avoid the computation of stiffness matrix due to the computational cost at the expense of improved convergence properties.

The iterative nature of a non-linear problem requires the use of tolerances to determine if the results are sufficiently converged.  For this, non-linear simulation cases need to specify the relative displacement and force residuals ::

    "NonlinearOptions" : {
        "DisplacementTolerance" : 1.0e-5,
        "ResidualTolerance" : 1.0e-4
    }

Methods to improve the properties of the Newton-Raphson could be implemented on top of the current non-linear solvers.


Non-linear implicit dynamic
===========================

The code to solve the implicit dynamic problem has been implemented but not yet tested.  If required, please provide a test case to the developers.
