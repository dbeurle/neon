
#pragma once

#include "femStaticMatrix.hpp"

#include "solver/time/GeneralisedTrapezoidal.hpp"

namespace neon::diffusion
{
class femDynamicMatrix : public femStaticMatrix
{
public:
    explicit femDynamicMatrix(fem_mesh& fem_mesh, json const& simulation_data);

    void solve() override final;

protected:
    /** Assembles the consistent (full) mass matrix */
    void assemble_mass();

protected:
    sparse_matrix M; //!< Consistent capacity matrix

    GeneralisedTrapezoidal time_solver;
};
}
