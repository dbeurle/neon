
#pragma once

#include "mesh/mechanical/solid/fem_mesh.hpp"
#include "numeric/sparse_matrix.hpp"
#include "solver/adaptive_time_step.hpp"

#include "io/FileIO.hpp"

namespace neon
{
class LinearSolver;
}

namespace neon::mechanical::solid
{
class fem_static_matrix
{
public:
    explicit fem_static_matrix(fem_mesh& mesh, json const& simulation);

    ~fem_static_matrix();

    void internal_restart(json const& solver_data, json const& new_increment_data);

    virtual void solve();

protected:
    /**
     * Compute the sparse pattern of the coefficient matrix but use the doublet
     * list as a means of forming the sparsity pattern.  This is a memory
     * intensive operation and should be replaced in the future
     */
    void compute_sparsity_pattern();

    void compute_internal_force();

    /** Gathers the external force contributions to the system of equations */
    void compute_external_force();

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
    void enforce_dirichlet_conditions(sparse_matrix& A, vector& b) const;

    /** Move the nodes on the mesh for the Dirichlet boundary */
    void apply_displacement_boundaries();

    /** Equilibrium iteration convergence criteria */
    bool is_iteration_converged() const;

    /** Pretty printer for the convergence of the Newton-Raphson solver */
    void print_convergence_progress() const;

    void update_relative_norms();

private:
    void perform_equilibrium_iterations();

protected:
    fem_mesh& mesh;

    FileIO<fem_mesh> io;

    adaptive_time_step adaptive_load;

    bool is_sparsity_computed = false;

    double residual_tolerance = 1.0e-3;
    double displacement_tolerance = 1.0e-3;

    int maximum_iterations = 10;

    double relative_displacement_norm;
    double relative_force_norm;

    sparse_matrix Kt; //!< Tangent matrix stiffness
    vector fint;      //!< Internal force vector
    vector fext;      //!< External force vector

    vector displacement;     //!< Displacement vector
    vector displacement_old; //!< Last displacement vector
    vector delta_d;          //!< Incremental displacement vector

    vector minus_residual; //!< Minus residual vector

    std::unique_ptr<LinearSolver> linear_solver;
};

inline bool fem_static_matrix::is_iteration_converged() const
{
    return relative_displacement_norm <= displacement_tolerance
           && relative_force_norm <= residual_tolerance;
}
}
