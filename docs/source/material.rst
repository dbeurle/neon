Material properties
===================

Constitutive models are used to model a material and require a number of material parameters.  These properties are tied closely to the constitutive model.  As an example, it is possible that multiple bodies have the same material parameters but different constitutive models.  To handle this particular case, each ``Material`` definition has a unique identifier ``Name``.  For example ::

    "Material" : [{
        "Name" : "steel"
    },
    {
        "Name" : "aluminium"
    }]

where the ``Name`` field acts as the unique key for each material.  An error will be produced for duplicated names.  This material field can have properties ::

    "Material" : [{
        "Name" : "steel",
        "ElasticModulus" : 134.0e3,
        "PoissonsRatio" : 0.3
    }]

and it is up to the constitutive model to check if the required material properties exist inside the material class.  This means you can include `extra` material properties and only the required values are used.  The user should check the constitutive model definition for the list of required material parameters.
