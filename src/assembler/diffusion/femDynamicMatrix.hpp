
#pragma once

#include "femStaticMatrix.hpp"

#include "solver/time/GeneralisedTrapezoidal.hpp"

namespace neon::diffusion
{
class femDynamicMatrix : public femStaticMatrix
{
public:
    explicit femDynamicMatrix(femMesh& fem_mesh,
                              Visualisation&& visualisation,
                              Json::Value const& solver_data,
                              Json::Value const& time_data);

    ~femDynamicMatrix() = default;

    void solve() override final;

protected:
    void assemble_mass();

protected:
    Vector M; //!< Diagonal mass matrix

    GeneralisedTrapezoidal trapezoidal;
};
}
