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
``"Traction"``  Specified with degree of freedom ``"x"``, ``"y"`` and/or ``"z"``
``"BodyForce"`` Specified with degree of freedom ``"x"``, ``"y"`` and/or ``"z"``
``"Pressure"``  Specified by ``"Value" : []``
=============== ============================================

Thermal
~~~~~~~

Thermal analysis specifies a heat flux which is defined to the normal of the surface.  In this case only a scalar value needs to be specified.

=============== ============================================
Keyword         Specification
=============== ============================================
``"HeatFlux"``  Specified by ``"Value" : []``
=============== ============================================

A special case occurs for the thermal analysis as a mixed (or Robin type) boundary condition.  This type is a combination of a Neumann and a Dirichlet condition and is used to model natural convection and is known as Newton Cooling.  This is specified by ::

    {
        "Name" : "fins",
        "Type" : "NewtonCooling",
        "Time" : [0.0, 5000.0],
        "HeatTransferCoefficient" : [50.0, 50.0],
        "AmbientTemperature" : [323.0, 323.0]
    }
