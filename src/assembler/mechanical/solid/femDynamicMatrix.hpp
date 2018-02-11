
#pragma once

#include "femStaticMatrix.hpp"

#include "solver/time/newmark_beta_integrator.hpp"

namespace neon::mechanical::solid
{
class femDynamicMatrix : public femStaticMatrix
{
public:
    explicit femDynamicMatrix(fem_mesh& fem_mesh, json const& simulation);

    ~femDynamicMatrix() = default;

    void solve() override final;

protected:
    void assemble_mass();

private:
    void perform_equilibrium_iterations();

protected:
    vector M; //!< Diagonal mass matrix

    vector a; //!< Nodal acceleration
    vector v; //!< Nodal velocity

    newmark_beta_integrator newmark;
};
}
