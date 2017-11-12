
#pragma once

#include "mesh/solid/femMesh.hpp"
#include "numeric/SparseTypes.hpp"
#include "solver/AdaptiveLoadStep.hpp"

#include "io/FileIO.hpp"

namespace neon
{
class LinearSolver;
}

namespace neon::mech::solid
{
class femStaticMatrix
{
public:
    explicit femStaticMatrix(femMesh& fem_mesh, Json::Value const& simulation);

    ~femStaticMatrix();

    void internal_restart(Json::Value const& solver_data, Json::Value const& new_increment_data);

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
     * Computes the external force contributions to the system of equations.
     * \sa NonFollowerLoadBoundary
     */
    void compute_external_force(double const time_step);

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
    bool is_iteration_converged() const;

    /** Pretty printer for the convergence of the Newton-Raphson solver */
    void print_convergence_progress() const;

    void update_relative_norms(Vector const& d_new, Vector const& delta_d, Vector const& residual);

private:
    void perform_equilibrium_iterations();

protected:
    femMesh& fem_mesh;

    FileIO file_io;

    AdaptiveLoadStep adaptive_load;

    bool is_sparsity_computed = false;

    double residual_tolerance = 1.0e-5;
    double displacement_tolerance = 1.0e-5;

    double relative_displacement_norm;
    double relative_force_norm;

    SparseMatrix Kt; //!< Tangent matrix stiffness
    Vector fint;     //!< Internal force vector
    Vector fext;     //!< External force vector
    Vector d;        //!< Displacement vector

    std::unique_ptr<LinearSolver> linear_solver;
};

inline bool femStaticMatrix::is_iteration_converged() const
{
    return relative_displacement_norm <= displacement_tolerance
           && relative_force_norm <= residual_tolerance;
}
}
