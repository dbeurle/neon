
#pragma once

#include "AbstractModule.hpp"

#include "assembler/mechanical/solid/femStaticMatrix.hpp"
#include "mesh/mechanical/solid/femMesh.hpp"

#include "assembler/mechanical/plane/femStaticMatrix.hpp"
#include "mesh/mechanical/plane/femMesh.hpp"

namespace neon
{
//! This namespace groups together all of the classes and functions associated
//! with solid mechanics finite elements.  The sub-namespaces implement mechanics
//! such as solid, plate, beam, axisymmetric etc
namespace mechanical
{
//! This namespace groups together all of the classes and functions associated
//! with three-dimensional solid mechanics finite elements.  These include
//! constitutive models, matrix systems, meshes, element stiffness matrices etc.
namespace solid
{
}
}

/**
 * SolidMechanicsModule is responsible for handling the setup and simulation of the class
 * of three dimensional solid mechanics problems.
 */
class SolidMechanicsModule : public AbstractModule
{
public:
    SolidMechanicsModule(BasicMesh const& mesh,
                         json const& material,
                         json const& simulation);

    virtual ~SolidMechanicsModule() = default;

    SolidMechanicsModule(SolidMechanicsModule const&) = delete;

    SolidMechanicsModule(SolidMechanicsModule&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    mechanical::solid::femMesh fem_mesh;           //!< Mesh with the solid routines
    mechanical::solid::femStaticMatrix fem_matrix; //!< Nonlinear solver routines
};

class PlaneMechanicsModule : public AbstractModule
{
public:
    PlaneMechanicsModule(BasicMesh const& mesh,
                         json const& material,
                         json const& simulation);

    virtual ~PlaneMechanicsModule() = default;

    PlaneMechanicsModule(PlaneMechanicsModule const&) = delete;

    PlaneMechanicsModule(PlaneMechanicsModule&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    mechanical::plane::femMesh fem_mesh;           //!< Mesh with the solid routines
    mechanical::plane::femStaticMatrix fem_matrix; //!< Nonlinear solver routines
};
}
