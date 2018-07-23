
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

/// beam_module holds the mesh and the matrix assembler for a linear C0
/// finite element beam formulation
class beam_module : public abstract_module
{
public:
    using mesh_type = mechanical::beam::fem_mesh;
    using matrix_type = fem::linear_static_matrix<mesh_type>;

public:
    explicit beam_module(basic_mesh const& mesh,
                         json const& material_data,
                         json const& profile_data,
                         json const& simulation_data);

    virtual ~beam_module() = default;

    beam_module(beam_module const&) = delete;

    beam_module(beam_module&&) = default;

    void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Beam finite element mesh
    mesh_type fem_mesh;
    /// Linear solver routines
    matrix_type fem_matrix;
};
}
