
#pragma once

#include "AbstractModule.hpp"

#include "assembler/mechanical/plane/femStaticMatrix.hpp"
#include "mesh/mechanical/plane/fem_mesh.hpp"

namespace neon
{
namespace mechanical
{
/// This namespace groups together all of the classes and functions associated
/// with plane strain finite elements.  These include constitutive models,
/// matrix systems, meshes, element stiffness matrices etc.
namespace plane_strain
{
}
}

class plane_strain_module : public AbstractModule
{
public:
    plane_strain_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~plane_strain_module() = default;

    plane_strain_module(plane_strain_module const&) = delete;

    plane_strain_module(plane_strain_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    mechanical::plane::fem_mesh fem_mesh;           //!< Mesh with the solid routines
    mechanical::plane::femStaticMatrix fem_matrix; //!< Nonlinear solver routines
};
}
