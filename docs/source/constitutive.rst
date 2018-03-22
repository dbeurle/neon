*******************
Constitutive models
*******************

Constitutive models play an important role in the solution of the equations of motion.  They are required to close the system such that the system is solvable.  These models relate quantities of displacement to the stress in the body.  Each model of material deformation requires some material properties and can effect the solution time, stability and accuracy of the overall finite element simulation.

All of the units are assumed to be consistent and it is the user's responsibility to make sure this happens.

Isotropic linear elasticity
===========================

The simplest form of constitutive model for linear materials is the isotropic linear elasticity model.  This is specified using the ``"IsotropicLinearElasticity"``.  The required parameters are

 * ``"PoissonsRatio"`` for the Poisson's Ratio of the material
 * ``"ElasticModulus"`` for the elasticity modulus (slope of stress-strain curve)

In addition, the elasticity constants can also be specified using the ``"ShearModulus"`` and the ``"BulkModulus"`` if this is more convenient to the user.

This model can be used by specifying ::

    "ConstitutiveModel" : {
        "Name" : "IsotropicLinearElasticity"
    }

Small strain J2 plasticity
==========================

The small strain J2 plasticity model introduces material non-linearity into the finite element formulation.  This model extends the ``"IsotropicLinearElasticity"`` model to include:

* ``"YieldStress"`` stress when the material begins plastic flow
* ``"IsotropicHardeningModulus"`` the hardening modulus of the material

This model is used by specifying ::

    "ConstitutiveModel" : {
        "Name" : "J2Plasticity",
        "FiniteStrain" : false
    }

where the finite strain version of the J2 model is not required.  Further modifications to the J2 model are given in subsequent subsections.

Isotropic Chaboche damage
~~~~~~~~~~~~~~~~~~~~~~~~~

TODO Fill out model details.

Neo-Hooke elasticity
====================

A common hyperelastic material model is the Neo-Hooke model which is a analogous to the linear elasticity model for finite strains.  The model is based on the compressible formulation is not incompressible.  Only two material parameters required are

    * ``"PoissonsRatio"`` for the Poisson's Ratio of the material
    * ``"ElasticModulus"`` for the elasticity modulus (slope of stress-strain curve)

alternatively the material parameters

    * ``"BulkModulus"``
    * ``"ShearModulus"``

can be used, where the ``"BulkModulus"`` acts as the penalty parameter for volumetric strain.

This model is used by specifying ::

    "ConstitutiveModel" : {
        "Name" : "Neohooke"
    }

Note that this constitutive model invokes the finite strain solver that requires an incremental approach and additional Newton-Raphson solver iterations.

Affine microsphere
==================

The affine microsphere model uses a Langevin statistical mechanics model of the polymer chains and derives a continuum model suitable for a finite element implementation.  This model uses numerical integration over a unit sphere to compute this homogenisation and requires a quadrature rule to be specified.  Note that higher order rules increase the computation time.  This model is suited to larger deformations than the Neo-Hooke model.

Material parameters are

    * ``"BulkModulus"``
    * ``"ShearModulus"``
    * ``"SegmentsPerChain"`` Number of segments per chain (experimental)

To specify the unit sphere quadrature rules

    * ``"BO21"`` 21 point symmetric scheme
    * ``"BO33"`` 33 point scheme
    * ``"FM900"`` 900 point scheme.  Note this is computationally expensive

The model is used by specifying ::

    "ConstitutiveModel" : {
        "Name" : "Microsphere",
        "Type" : "Affine",
        "Statistics" : "Langevin",
        "Quadrature" : "BO21"
    }


Gaussian affine microsphere
===========================

The Gaussian affine microsphere model re-derives the affine microsphere model using a Gaussian chain description.  This significantly reduces complexity of the model.

Material parameters are

    * ``"BulkModulus"``
    * ``"ShearModulus"``
    * ``"SegmentsPerChain"`` Number of segments per chain (not required)

To specify the unit sphere quadrature rules

    * ``"BO21"`` 21 point symmetric scheme
    * ``"BO33"`` 33 point scheme
    * ``"FM900"`` 900 point scheme.  Note this is computationally expensive

The model is used by specifying ::

    "ConstitutiveModel" : {
        "Name" : "Microsphere",
        "Type" : "Affine",
        "Statistics" : "Gaussian",
        "Quadrature" : "BO21"
    }
