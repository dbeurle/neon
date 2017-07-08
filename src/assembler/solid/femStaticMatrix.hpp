
#pragma once

#include "mesh/solid/femMesh.hpp"
#include "numeric/SparseTypes.hpp"
#include "solver/AdaptiveIncrement.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearSolver;
}

namespace neon::solid
{
class femStaticMatrix
{
public:
    femStaticMatrix(femMesh& fem_mesh, Json::Value const& solver_data);

    ~femStaticMatrix();

    virtual void solve();

protected:
    /**
     * Compute the sparse pattern of the coefficient matrix but use the doublet
     * list as a means of forming the sparsity pattern.  This is a memory
     * intensive operation and should be replaced in the future
     */
    void compute_sparsity_pattern();

    void compute_internal_force();

    /**
     * Assembles the material and geometric matrices, checking for allocation
     * already performed
     */
    void assemble_stiffness();

    /**
     * Apply dirichlet conditions to the system defined by A, x, and b.
     * This method sets the incremental displacements to zero for the given
     * load increment such that incremental displacements are zero
     */
    void enforce_dirichlet_conditions(SparseMatrix& A, Vector& x, Vector& b);

    void apply_displacement_boundaries(double const load_factor = 1.0);

private:
    void perform_equilibrium_iterations();

protected:
    femMesh& fem_mesh;

    AdaptiveIncrement adaptive_load;

    bool is_sparsity_computed = false;

    SparseMatrix Kt; //!< Tangent matrix stiffness
    Vector fint;     //!< Internal force vector
    Vector d;        //!< Displacement

    std::unique_ptr<LinearSolver> linear_solver;
};
}
