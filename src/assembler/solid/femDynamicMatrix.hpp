
#pragma once

#include "femStaticMatrix.hpp"

#include "solver/time/NewmarkBeta.hpp"

namespace neon::solid
{
class femDynamicMatrix : public femStaticMatrix
{
public:
    explicit femDynamicMatrix(femMesh& fem_mesh,
                              Visualisation&& visualisation,
                              Json::Value const& solver_data,
                              Json::Value const& nonlinear_data,
                              Json::Value const& time_data);

    ~femDynamicMatrix() = default;

    void solve() override final;

protected:
    void assemble_mass();

private:
    void perform_equilibrium_iterations();

protected:
    Vector M; //!< Diagonal mass matrix

    Vector a; //!< Nodal acceleration
    Vector v; //!< Nodal velocity

    NewmarkBeta newmark;
};
}
