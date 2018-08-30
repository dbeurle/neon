*******************
Boundary Conditions
*******************

An essential part of defining the model is the specification of appropriate boundary conditions.  Boundary conditions are applied to the mesh and included in the mathematical description of the system.  There are a few different types of boundary conditions, which are specialised for each particular type of simulation.

Each boundary condition is applied over a specified time.  That means, for each boundary condition the user must specify *at least* two data points inside a vector type ::

    "boundaries" : [{
        "name" : "x-sym",
        "type" : "displacement",
        "time" : [0.0, 1.0],
        "x" : [0.0, 0.0]
    }]

where multiple boundary conditions can be added by inserting extra elements into the array.

The specified time and values are automatically linearly interpolated if a time step happens to fall between two values.  For convenience the adaptive time step algorithm will always produce a solution at the specified time points to enable easy experimental comparisons where only some discrete time points are evaluated.

Since boundary conditions can vary over time, it is possible to specify cyclic type boundary conditions that follow a periodic curve.  Currently, only a sinusoidal periodic type is implemented and it can be specified with ::

    "boundaries" : [{
        "name" : "x-sym",
        "type" : "displacement",
        "generate_type" : "sinusoidal",
        "x" : [0.05, 0.10],
        "period" : [10,10],
        "phase" : [0,0],
        "number_of_cycles" : [2,1]
    }]


This may be used to generate block loading where each block follows a specific sinusoidal wave that corresponds to one entry of ``"x"``, ``"period"``, ``"phase"`` and ``"number_of_cycles"``. In this case the time step is taken to be equal to ::

    "time" : {
        "period" : 30,
        "increments" : {
            "initial" : 0.1
        }

and the time steps are generated automatically based on the given boundary condition parameters.


Essential (Dirichlet) Type
==========================

These boundary conditions are required to ensure there is no rigid body motion of the system under load.

=============== ============================================
Simulation type Specification
=============== ============================================
Solid mechanics ``"displacement"`` specified with degree of freedom ``"x"``, ``"y"`` or ``"z"``
Thermal         ``"temperature"`` specified by ``"value" : []``
=============== ============================================

Natural (Neumann) Type
======================

In contrast to Dirichlet boundary conditions, Neumann boundary conditions specify a derivative over a boundary.  These can be used to model tractions, pressure, heat flux, gravitational loads and internal heat generation.  Because of their diverse meaning, this section will be split into each simulation type possible.

Solid Mechanics
~~~~~~~~~~~~~~~

================  ============================================
Keyword           Specification
================  ============================================
``"traction"``    Specified with degree of freedom ``"x"``, ``"y"`` and/or ``"z"``
``"body_force"``  Specified with degree of freedom ``"x"``, ``"y"`` and/or ``"z"``
``"pressure"``    Specified by ``"value" : []``
================  ============================================

Note here, that ``"pressure"`` is a non-follower load and does not rotate with the deformation of the body.

Thermal
~~~~~~~

Thermal analysis specifies a heat flux which is defined to the normal of the surface.  In this case only a scalar value needs to be specified.

=============== ============================================
Keyword         Specification
=============== ============================================
``"flux"``      Specified by ``"value" : []``
=============== ============================================

A special case occurs for the thermal analysis as a mixed (or Robin type) boundary condition.  This type is a combination of a Neumann and a Dirichlet condition and is used to model natural convection and is known as Newton Cooling.  This is specified by ::

    {
        "name" : "fins",
        "type" : "newton_cooling",
        "time" : [0.0, 5000.0],
        "heat_transfer_coefficient" : [50.0, 50.0],
        "ambient_temperature" : [323.0, 323.0]
    }
