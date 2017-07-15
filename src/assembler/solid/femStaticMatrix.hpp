
#pragma once

#include "mesh/solid/femMesh.hpp"
#include "numeric/SparseTypes.hpp"
#include "solver/AdaptiveLoadStep.hpp"
#include "visualisation/Visualisation.hpp"

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
    explicit femStaticMatrix(femMesh& fem_mesh,
                             Visualisation&& visualisation,
                             Json::Value const& solver_data,
                             Json::Value const& increment_data);

    ~femStaticMatrix();

    void internal_restart(Json::Value const& new_increment_data);

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

    /** Move the nodes on the mesh for the Dirichlet boundary */
    void apply_displacement_boundaries();

    /** Equilibrium iteration convergence criteria */
    bool is_converged(double const inc_disp_norm, double const residual_norm) const;

    /** Pretty printer for the convergence of the Newton-Raphson solver */
    void print_convergence_progress(double const delta_d_norm,
                                    double const residual_norm) const;

private:
    void perform_equilibrium_iterations();

protected:
    femMesh& fem_mesh;

    Visualisation visualisation;

    AdaptiveLoadStep adaptive_load;

    bool is_sparsity_computed = false;

    double residual_tolerance = 1.0e-4;
    double displacement_tolerance = 1.0e-6;

    SparseMatrix Kt; //!< Tangent matrix stiffness
    Vector fint;     //!< Internal force vector
    Vector d;        //!< Displacement

    std::unique_ptr<LinearSolver> linear_solver;
};

inline bool femStaticMatrix::is_converged(double const inc_disp_norm,
                                          double const residual_norm) const
{
    return inc_disp_norm <= displacement_tolerance && residual_norm <= residual_tolerance;
}
}
