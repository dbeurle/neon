
#pragma once

#include "abstract_module.hpp"

#include "assembler/linear_static_matrix.hpp"
#include "mesh/mechanical/beam/fem_mesh.hpp"

namespace neon
{
namespace mechanical
{
/// This namespace groups together all of the classes and functions associated
/// with beam formulation in finite elements.  These include constitutive models,
/// matrix systems, meshes, element stiffness matrices etc.
namespace beam
{
}
}

class beam_module : public abstract_module
{
public:
    using mesh_type = mechanical::beam::fem_mesh;
    using matrix_type = fem::linear_static_matrix<mesh_type>;

public:
    beam_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~beam_module() = default;

    beam_module(beam_module const&) = delete;

    beam_module(beam_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Mesh for plane types
    mesh_type fem_mesh;
    /// Nonlinear solver routines
    matrix_type fem_matrix;
};
}
