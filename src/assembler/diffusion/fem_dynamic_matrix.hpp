
#pragma once

#include "fem_static_matrix.hpp"

#include "solver/time/trapezoidal_integrator.hpp"

namespace neon::diffusion
{
class fem_dynamic_matrix : public fem_static_matrix
{
public:
    explicit fem_dynamic_matrix(fem_mesh& fem_mesh, json const& simulation_data);

    void solve() override final;

protected:
    /** Assembles the consistent (full) mass matrix */
    void assemble_mass();

protected:
    sparse_matrix M; //!< Consistent capacity matrix

    trapezoidal_integrator time_solver;
};
}
