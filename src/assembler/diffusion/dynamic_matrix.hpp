
#pragma once

#include "static_matrix.hpp"

#include "solver/time/trapezoidal_integrator.hpp"

namespace neon::diffusion
{
class dynamic_matrix : public static_matrix
{
public:
    using static_matrix::mesh_type;

public:
    explicit dynamic_matrix(mesh_type& mesh, json const& simulation_data);

    void solve() override final;

protected:
    /// Assembles the consistent (full) mass matrix
    void assemble_mass();

protected:
    /// Consistent capacity matrix
    sparse_matrix M;

    trapezoidal_integrator time_solver;
};
}
