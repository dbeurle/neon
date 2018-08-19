*******************
Constitutive Models
*******************

Constitutive models play an important role in the solution of the equations of motion.  They are required to close the system such that the system is solvable.  These models relate quantities of displacement to the stress in the body.  Each model of material deformation requires some material properties and can effect the solution time, stability and accuracy of the overall finite element simulation.

All of the units are assumed to be consistent and it is the user's responsibility to make sure this happens.

Isotropic Linear Elasticity
===========================

The simplest form of constitutive model for linear materials is the isotropic linear elasticity model.  This is specified using  ``"IsotropicLinearElasticity"``.  The required parameters are

 * ``"PoissonsRatio"`` for the Poisson's Ratio of the material
 * ``"ElasticModulus"`` for the elasticity modulus (slope of stress-strain curve)

In addition, the elasticity constants can also be specified using  ``"ShearModulus"`` and ``"BulkModulus"`` if this is more convenient for the user.

The model is specified by ::

    "ConstitutiveModel" : {
        "Name" : "IsotropicLinearElasticity"
    }

Small Strain J2 Plasticity
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

Isotropic Chaboche Damage
~~~~~~~~~~~~~~~~~~~~~~~~~

This non-linear isotropic ductile damage model is based on the ``"J2Plasticity"`` model. It is mainly controlled by the following material parameters:

* ``"KinematicHardeningModulus"`` the kinematic hardening modulus of the material
* ``"SofteningMultiplier"``  a scaling factor multiplied by the back stress in the plastic yield function
* ``"PlasticityViscousExponent"`` :math:`n`: relates the plastic multiplier to the plastic yield function exponentially
* ``"PlasticityViscousDenominator"`` :math:`K`: scales the plastic yield function :math:`\dot{\lambda}_{vp} = \langle \frac{f_{vp}}{K} \rangle^{n}`
* ``"DamageViscousExponent"`` :math:`s`: relates the damage multiplier to the damage yield function exponentially
* ``"DamageViscousDenominator"`` :math:`S`: scales the damage yield function :math:`\dot{\lambda}_{d} = \langle \frac{f_{d}}{S} \rangle^{s}`

This model is used by specifying ::

    "ConstitutiveModel" : {
          "Name" : "J2Plasticity",
          "Damage" : "IsotropicChaboche",
          "FiniteStrain" : false
      }

This constitutive model is implemented only in an infinitesimal strain framework.

Neo-Hooke Elasticity
====================

The Neo-Hooke model is a common hyperelastic material model which can be used to model non-linear elasticity behaviour.  Only two material parameters required, which are

    * ``"PoissonsRatio"`` for the Poisson's ratio of the material
    * ``"ElasticModulus"`` for the elasticity modulus (slope of stress-strain curve)

alternatively the material parameters

    * ``"BulkModulus"``
    * ``"ShearModulus"``

can be used, where the ``"BulkModulus"`` acts as the penalty parameter for volumetric strain.

This model is used by specifying ::

    "ConstitutiveModel" : {
        "Name" : "Neohooke"
    }

Note that this constitutive model invokes the finite strain solver, which requires an incremental approach and additional Newton-Raphson solver iterations.

Affine Microsphere
==================

The affine microsphere model uses a Langevin statistical mechanics model for polymer chains and derives a continuum model suitable for a finite element implementation.  This model uses numerical integration over a unit sphere to compute this homogenisation and requires a quadrature rule to be specified.  Note that higher order rules increase the computation time.  This model is suited to larger deformations than the Neo-Hooke model.

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


Gaussian Affine Microsphere
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
