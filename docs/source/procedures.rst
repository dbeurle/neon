Procedures
==========

There are several procedures availble to handle structural analysis in the finite element framework.  The most general procedure is the continuum theory, which solves the momentum equation over a three dimensional body.  However, the computational cost of this can be high which led to the development of specialised beam and shell theory.  These reduced models alleviate the burden of a three-dimensional simulation with some restrictions.

In this section, an outline of specialised theory in ``neon`` is documented along with the modelling assumption the theory contains.

Beam
----

When approximating the structural member with beam theory, additional information regarding the cross-sectional area and the second moment of area must be given in order find an approximate solution.  This is accomplished for all beam theories by specifying a ``profile`` and the orientation of this profile in space for every finite element.

Each element group with the same profile is grouped together inside a ``mesh`` object in the input file ::

    "mesh" : [{
        ...
        "section" : [{
            "name" : "<name of element group in mesh file>",
            "profile" : "<name of profile definition>",
            "tangent" : [0.0, 1.0, 0.0],
            "normal" : [1.0, 0.0, 0.0]
        },
        ...
        ],
    }]

Linear (C0)
~~~~~~~~~~~

In order to avoid the requirement of higher order shape functions (see Beam Theory (C1)), a reduced form of the momentum equation is solved rather than discretising the Euler-Bernouilli fourth order differential equation.  This involves the computation of four separate stiffness matrices; axial, torsional, bending and shear.  A specialised integration scheme is automatically applied to this analysis to avoid shear locking through under integration of selected stiffness matrices.



Euler-Bernouilli (C1)
~~~~~~~~~~~~~~~~~~~~~

If the Euler-Bernioulli beam theory is directly formulated in the finite element method, the resulting shape functions are required to be C1 continuous, resulting in a higher order interpolation.  This theory will be implemented after the C0 theory.

Continuum
---------

The evaluation of the integrals in the finite element method use numerical quadrature rules specialised for each element.  For computational efficiency reasons, these could be under-integrated or fully-integrated.  However, in the current state reduced integration rules produce a rank-deficient element stiffness matrix and are therefore not recommended for use until stabilisation is applied.  For each mesh type, the element options can be specified ::

    "ElementOptions" : {
        "Quadrature" : "Full"
    }

with reduced integration selected when ``"Quadrature" : "Reduced"``.
