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

where the finite strain version of the J2 model is not required.
