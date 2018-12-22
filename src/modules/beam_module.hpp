
#pragma once

/// @file

#include "abstract_module.hpp"

#include "assembler/linear_static_matrix.hpp"
#include "mesh/mechanics/beam/mesh.hpp"

namespace neon
{
namespace mechanics
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
    using mesh_type = mechanics::beam::mesh;
    using matrix_type = linear_static_matrix<mesh_type>;

public:
    explicit beam_module(basic_mesh const& mesh,
                         json const& material_data,
                         json const& simulation_data,
                         std::map<std::string, std::unique_ptr<geometry::profile>> const& profile_store);

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
