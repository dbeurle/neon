
#pragma once

#include "mesh/solid/femMesh.hpp"
#include "numeric/SparseTypes.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearSolver;
}

namespace neon::solid
{
class femMatrix
{
public:
    femMatrix(femMesh& fem_mesh, Json::Value const& solver_data);

    ~femMatrix();

    void solve();

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

    void apply_displacements(double const load_factor = 1.0);

protected:
    femMesh& fem_mesh;

    bool is_sparsity_computed = false;

    SparseMatrix Kt; //!< Tangent matrix stiffness
    Vector fint;     //!< Internal force vector

    std::unique_ptr<LinearSolver> linear_solver;
};
}
