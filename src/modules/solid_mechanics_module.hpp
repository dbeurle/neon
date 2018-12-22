
#pragma once

/// @file

#include "abstract_module.hpp"

#include "assembler/mechanics/static_matrix.hpp"
#include "assembler/mechanics/latin_matrix.hpp"

#include "assembler/mechanics/linear_buckling_matrix.hpp"
#include "assembler/mechanics/natural_frequency_matrix.hpp"

#include "mesh/mechanics/solid/mesh.hpp"

namespace neon
{
/// This namespace groups together all of the classes and functions associated
/// with solid mechanics finite elements.  The sub-namespaces implement mechanics
/// such as solid, plate, beam, axisymmetric etc
namespace mechanics
{
/// This namespace groups together all of the classes and functions associated
/// with three-dimensional solid mechanics finite elements.  These include
/// constitutive models, matrix systems, meshes, element stiffness matrices etc.
namespace solid
{
}
}

/// solid_mechanics_module is responsible for handling the setup and simulation
///  of the class for three dimensional solid mechanics problems.
template <typename matrix_type>
class solid_mechanics_module : public abstract_module
{
public:
    using mesh_type = typename matrix_type::mesh_type;

public:
    solid_mechanics_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~solid_mechanics_module() = default;

    solid_mechanics_module(solid_mechanics_module const&) = delete;

    solid_mechanics_module(solid_mechanics_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Mesh for solid types
    mesh_type fem_mesh;
    /// Nonlinear solver routines
    matrix_type fem_matrix;
};
extern template class solid_mechanics_module<
    mechanics::static_matrix<mechanics::solid::mesh<mechanics::solid::submesh>>>;
extern template class solid_mechanics_module<
    mechanics::latin_matrix<mechanics::solid::mesh<mechanics::solid::latin_submesh>>>;

/// linear_buckling_module is responsible for handling the setup
/// and simulation of the class for three dimensional linear (eigenvalue)
/// solid mechanics buckling problems.
class linear_buckling_module : public abstract_module
{
public:
    using mesh_type = mechanics::solid::mesh<mechanics::solid::submesh>;
    using matrix_type = mechanics::linear_buckling_matrix<mesh_type>;

public:
    linear_buckling_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~linear_buckling_module() = default;

    linear_buckling_module(linear_buckling_module const&) = delete;

    linear_buckling_module(linear_buckling_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Mesh with the solid routines
    mesh_type fem_mesh;
    /// Linear eigenvalue buckling solver
    matrix_type fem_matrix;
};

/// natural_frequency_module is responsible for handling the setup
/// and simulation of the class for eigenvalue solid mechanics natural frequency
/// problems
class natural_frequency_module : public abstract_module
{
public:
    using mesh_type = mechanics::solid::mesh<mechanics::solid::submesh>;
    using matrix_type = mechanics::natural_frequency_matrix<mesh_type>;

public:
    natural_frequency_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~natural_frequency_module() = default;

    natural_frequency_module(natural_frequency_module const&) = delete;

    natural_frequency_module(natural_frequency_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Linear eigenvalue buckling solver
    matrix_type fem_matrix;
};
}
