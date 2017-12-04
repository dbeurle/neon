
#pragma once

#include "AbstractModule.hpp"

#include "assembler/solid/femStaticMatrix.hpp"
#include "mesh/mech/solid/femMesh.hpp"

namespace neon
{
//! This namespace groups together all of the classes and functions associated
//! with solid mechanics finite elements.  The sub-namespaces implement mechanics
//! such as solid, plate, beam, axisymmetric etc
namespace mech
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
                         Json::Value const& material,
                         Json::Value const& simulation);

    virtual ~SolidMechanicsModule() = default;

    SolidMechanicsModule(SolidMechanicsModule const&) = delete;

    SolidMechanicsModule(SolidMechanicsModule&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    mech::solid::femMesh fem_mesh;           //!< Mesh with the solid routines
    mech::solid::femStaticMatrix fem_matrix; //!< Nonlinear solver routines
};
}
