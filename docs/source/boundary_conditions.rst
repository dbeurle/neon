*******************
Boundary conditions
*******************

An essential part of defining the model is the specification of appropriate boundary conditions.  Boundary conditions are applied to the mesh and included in the mathematical description of the system.  There are a few different types of boundary conditions, which are specialised for each particular type of simulation.

Each boundary condition is applied over a specified time.  That means, for each boundary condition the user must specify *at least* two data points inside a vector type ::

    "BoundaryConditions" : [{
        "Name" : "x-sym",
        "Type" : "Displacement",
        "Time" : [0.0, 1.0],
        "x" : [0.0, 0.0]
    }]

where multiple boundary conditions can be added by inserting extra elements into the array.

The specified time and values are automatically linearly interpolated if a time step happens to fall between two values.  For convenience the adaptive time step algorithm will always produce a solution at the specified time points to easy experimental comparisons.

Since boundary conditions can vary over time, it is possible to specify cyclic type boundary conditions that follow a periodic curve.  This can be specified with ::

    "Periodic"

Essential (Dirichlet) type
==========================

These boundary conditions are required to ensure there is no rigid body motion of the system under load.

=============== ============================================
Simulation type Specification
=============== ============================================
Solid mechanics ``"Displacement"`` specified with degree of freedom ``"x"``, ``"y"`` or ``"z"``
Thermal         ``"Temperature"`` specified by ``"Value" : []``
=============== ============================================

Natural (Neumann) type
======================

In contrast to Dirichlet boundary conditions, Neumann boundary conditions specify a derivative over a boundary.  These can be used to model tractions, pressure, heat flux, gravitational loads and internal heat generation.  Because of their diverse meaning, this section will be split into each simulation type possible.

Solid mechanics
~~~~~~~~~~~~~~~

=============== ============================================
Keyword         Specification
=============== ============================================
``"Traction"``  Specified with degree of freedom ``"x"``, ``"y"`` or ``"z"``
``"BodyForce"`` Specified with degree of freedom ``"x"``, ``"y"`` or ``"z"``
``"Pressure"``  Specified by ``"Value" : []``
=============== ============================================
