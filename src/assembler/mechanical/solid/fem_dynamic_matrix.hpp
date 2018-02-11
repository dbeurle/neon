
#pragma once

#include "fem_static_matrix.hpp"

#include "solver/time/newmark_beta_integrator.hpp"

namespace neon::mechanical::solid
{
class fem_dynamic_matrix : public fem_static_matrix
{
public:
    explicit fem_dynamic_matrix(fem_mesh& fem_mesh, json const& simulation);

    ~fem_dynamic_matrix() = default;

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
