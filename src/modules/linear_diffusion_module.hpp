
#pragma once

#include "abstract_module.hpp"

#include "assembler/diffusion/fem_dynamic_matrix.hpp"
#include "assembler/diffusion/fem_static_matrix.hpp"
#include "mesh/diffusion/fem_mesh.hpp"

namespace neon
{
//! This namespace groups together all of the classes and functions associated
//! with three-dimensional diffusion type finite elements.  These include
//! constitutive models, matrix systems, meshes, element stiffness matrices etc.
namespace diffusion
{
}

/**
 * linear_diffusionModule is responsible for the construction and solution method
 * of a linear diffusion problem
 */
template <typename femMatrix_Tp>
class linear_diffusion_module : public abstract_module
{
public:
    explicit linear_diffusion_module(basic_mesh const& mesh,
                                     json const& material,
                                     json const& simulation)
        : fem_mesh(mesh, material, simulation["Mesh"][0]), fem_matrix(fem_mesh, simulation)
    {
    }

    virtual ~linear_diffusion_module() = default;

    void perform_simulation() override final { fem_matrix.solve(); }

protected:
    diffusion::fem_mesh fem_mesh;
    femMatrix_Tp fem_matrix;
};
}
