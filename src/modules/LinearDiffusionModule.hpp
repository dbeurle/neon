
#pragma once

#include "AbstractModule.hpp"

#include "assembler/diffusion/femStaticMatrix.hpp"
#include "mesh/diffusion/femMesh.hpp"

namespace neon
{
/**
 * LinearDiffusionModule is responsible for the construction and solution method
 * of a linear diffusion problem
 */
class LinearDiffusionModule : public AbstractModule
{
public:
    LinearDiffusionModule(BasicMesh const& mesh,
                          Json::Value const& material,
                          Json::Value const& simulation);

    void perform_simulation() override final { fem_matrix.solve(); }

protected:
    diffusion::femMesh fem_mesh;
    diffusion::femStaticMatrix fem_matrix;
};
}
