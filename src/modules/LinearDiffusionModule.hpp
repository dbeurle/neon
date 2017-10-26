
#pragma once

#include "AbstractModule.hpp"

#include "assembler/diffusion/femDynamicMatrix.hpp"
#include "assembler/diffusion/femStaticMatrix.hpp"
#include "mesh/diffusion/femMesh.hpp"

namespace neon
{
//! This namespace groups together all of the classes and functions associated
//! with three-dimensional diffusion type finite elements.  These include
//! constitutive models, matrix systems, meshes, element stiffness matrices etc.
namespace diffusion
{
}

/**
 * LinearDiffusionModule is responsible for the construction and solution method
 * of a linear diffusion problem
 */
template <typename femMatrix_Tp>
class LinearDiffusionModule : public AbstractModule
{
public:
    explicit LinearDiffusionModule(BasicMesh const& mesh,
                                   Json::Value const& material,
                                   Json::Value const& simulation)
        : fem_mesh(mesh, material, simulation["Mesh"][0]), fem_matrix(fem_mesh, simulation)
    {
    }

    virtual ~LinearDiffusionModule() = default;

    void perform_simulation() override final { fem_matrix.solve(); }

protected:
    diffusion::femMesh fem_mesh;
    femMatrix_Tp fem_matrix;
};
}
