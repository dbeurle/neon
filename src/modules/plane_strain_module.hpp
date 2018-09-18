
#pragma once

#include "abstract_module.hpp"

#include "assembler/mechanics/static_matrix.hpp"
#include "mesh/mechanics/plane/mesh.hpp"

namespace neon
{
namespace mechanics
{
/// This namespace groups together all of the classes and functions associated
/// with plane strain finite elements.  These include constitutive models,
/// matrix systems, meshes, element stiffness matrices etc.
namespace plane_strain
{
}
}

class plane_strain_module : public abstract_module
{
public:
    using mesh_type = mechanics::plane::mesh;
    using matrix_type = mechanics::static_matrix<mesh_type>;

public:
    plane_strain_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~plane_strain_module() = default;

    plane_strain_module(plane_strain_module const&) = delete;

    plane_strain_module(plane_strain_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Mesh for plane types
    mesh_type fem_mesh;
    /// Nonlinear solver routines
    matrix_type fem_matrix;
};
}
