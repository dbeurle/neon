Material Properties
===================

Constitutive models are used to model a material and require a number of material parameters.  These properties are tied closely to the constitutive model.  As an example, it is possible for multiple bodies to have the same material parameters but different constitutive models.  To handle this particular case, each ``materials`` definition has a unique identifier ``name``.  For example ::

    "materials" : [
    {
        "name" : "steel"
    },
    {
        "name" : "aluminium"
    }
    ]

where the ``name`` field acts as the unique key for each material.  An error will be produced for duplicated names.  This material field can have properties ::

    "materials" : [{
        "name" : "steel",
        "elastic_modulus" : 134.0e3,
        "poissons_ratio" : 0.3
    }]

and it is up to the constitutive model and the material model to check if the required material properties exist inside the material class.  This means you can include irrelevant material properties and only the required values are used.  The user should check the constitutive model definition for the list of required material parameters.
