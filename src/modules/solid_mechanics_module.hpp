
#pragma once

#include "abstract_module.hpp"

#include "assembler/mechanics/static_matrix.hpp"
#include "assembler/mechanics/buckling_matrix.hpp"

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
class solid_mechanics_module : public abstract_module
{
public:
    using mesh_type = mechanics::solid::mesh;
    using matrix_type = mechanics::fem_static_matrix<mesh_type>;

public:
    solid_mechanics_module(basic_mesh const& mesh, json const& material, json const& simulation);

    virtual ~solid_mechanics_module() = default;

    solid_mechanics_module(solid_mechanics_module const&) = delete;

    solid_mechanics_module(solid_mechanics_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Mesh for solid types
    mechanics::solid::mesh fem_mesh;
    /// Nonlinear solver routines
    matrix_type fem_matrix;
};

/// solid_mechanics_linear_buckling_module is responsible for handling the setup
/// and simulation of the class for three dimensional linear (eigenvalue)
/// solid mechanics buckling problems.
class solid_mechanics_linear_buckling_module : public abstract_module
{
public:
    using mesh_type = mechanics::solid::mesh;
    using matrix_type = mechanics::fem_buckling_matrix<mesh_type>;

public:
    solid_mechanics_linear_buckling_module(basic_mesh const& mesh,
                                           json const& material,
                                           json const& simulation);

    virtual ~solid_mechanics_linear_buckling_module() = default;

    solid_mechanics_linear_buckling_module(solid_mechanics_linear_buckling_module const&) = delete;

    solid_mechanics_linear_buckling_module(solid_mechanics_linear_buckling_module&&) = default;

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    /// Mesh with the solid routines
    mesh_type fem_mesh;
    /// Linear eigenvalue buckling solver
    matrix_type fem_matrix;
};
}
