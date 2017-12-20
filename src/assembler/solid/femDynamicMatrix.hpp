
#pragma once

#include "femStaticMatrix.hpp"

#include "solver/time/NewmarkBeta.hpp"

namespace neon::mechanical::solid
{
class femDynamicMatrix : public femStaticMatrix
{
public:
    explicit femDynamicMatrix(femMesh& fem_mesh, Json::Value const& simulation);

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

    NewmarkBeta newmark;
};
}
