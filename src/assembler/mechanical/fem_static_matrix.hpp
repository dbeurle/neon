
#pragma once

#include "assembler/sparsity_pattern.hpp"
#include "numeric/float_compare.hpp"
#include "exceptions.hpp"
#include "numeric/sparse_matrix.hpp"
#include "solver/adaptive_time_step.hpp"
#include "solver/linear/linear_solver_factory.hpp"
#include "io/file_output.hpp"
#include "io/json.hpp"

#include <chrono>
#include <memory>
#include <string>
#include <iostream>

#include <termcolor/termcolor.hpp>
#include <tbb/parallel_for.h>

namespace neon::mechanical
{
/// Generic static matrix designed for solid mechanics problems using the
/// Newton-Raphson method for the solution of nonlinear equations.
/// This class is responsible for the assembly of the process stiffness matrix,
/// residual force vector and solution of the incremental displacement.
template <class fem_mesh_type>
class fem_static_matrix
{
public:
    using mesh_type = fem_mesh_type;

public:
    explicit fem_static_matrix(mesh_type& fem_mesh, json const& simulation);

    /// Solve the nonlinear system of equations
    void solve();

protected:
    /// Gathers the internal force vector using the Cauchy stress
    void compute_internal_force();

    /// Gathers the external force contributions to the system of equations
    void compute_external_force();

    /// Assembles the material and geometric matrices, checking for allocation
    /// already performed
    void assemble_stiffness();

    /// Apply dirichlet conditions to the system defined by A, x, and b.
    /// This method sets the incremental displacements to zero for the given
    /// load increment such that incremental displacements are zero
    void enforce_dirichlet_conditions(sparse_matrix& A, vector& b) const;

    /// Move the nodes on the mesh for the Dirichlet boundary
    void apply_displacement_boundaries();

    /// Equilibrium iteration convergence criteria
    bool is_iteration_converged() const;

    /// Pretty printer for the convergence of the Newton-Raphson solver
    void print_convergence_progress() const;

    void update_relative_norms();

private:
    void perform_equilibrium_iterations();

protected:
    mesh_type& fem_mesh;

    file_output<mesh_type> io;

    adaptive_time_step adaptive_load;

    /// Cache the sparsity pattern
    bool is_sparsity_computed{false};

    double residual_tolerance{1.0e-3};
    double displacement_tolerance{1.0e-3};

    double relative_displacement_norm;
    double relative_force_norm;

    /// Tangent sparse stiffness matrix
    sparse_matrix Kt;
    /// Internal force vector
    vector f_int;
    /// External force vector
    vector f_ext;

    /// Displacement vector
    vector displacement;
    /// Last displacement vector
    vector displacement_old;
    /// Incremental displacement vector
    vector delta_d;

    /// Minus residual vector
    vector minus_residual;

    std::unique_ptr<linear_solver> solver;
};

template <class fem_mesh_type>
fem_static_matrix<fem_mesh_type>::fem_static_matrix(mesh_type& fem_mesh, json const& simulation)
    : fem_mesh(fem_mesh),
      io(simulation["Name"], simulation["Visualisation"], fem_mesh),
      adaptive_load(simulation["Time"], fem_mesh.time_history()),
      solver(make_linear_solver(simulation["LinearSolver"], fem_mesh.is_symmetric()))
{
    if (!simulation["NonlinearOptions"].count("DisplacementTolerance"))
    {
        throw std::domain_error("DisplacementTolerance not specified in "
                                "NonlinearOptions");
    }
    if (!simulation["NonlinearOptions"].count("ResidualTolerance"))
    {
        throw std::domain_error("ResidualTolerance not specified in "
                                "NonlinearOptions");
    }
    residual_tolerance = simulation["NonlinearOptions"]["ResidualTolerance"];
    displacement_tolerance = simulation["NonlinearOptions"]["DisplacementTolerance"];

    f_int = f_ext = displacement = displacement_old = delta_d = vector::Zero(fem_mesh.active_dofs());

    // Perform Newton-Raphson iterations
    std::cout << "\n"
              << std::string(4, ' ') << "Non-linear equation system has " << fem_mesh.active_dofs()
              << " degrees of freedom\n";
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::solve()
{
    try
    {
        // Initialise the mesh with zero displacements
        fem_mesh.update_internal_variables(displacement);

        while (!adaptive_load.is_fully_applied())
        {
            std::cout << "\n"
                      << std::string(4, ' ') << termcolor::magenta << termcolor::bold
                      << "Performing equilibrium iterations for time " << adaptive_load.step_time()
                      << termcolor::reset << std::endl;

            compute_external_force();

            perform_equilibrium_iterations();
        }
    }
    catch (computational_error& comp_error)
    {
        // A numerical error has been reported that is able to be recovered
        // by resetting the state
        std::cout << std::endl
                  << std::string(6, ' ') << termcolor::bold << termcolor::yellow
                  << comp_error.what() << termcolor::reset << std::endl;

        adaptive_load.update_convergence_state(false);
        fem_mesh.save_internal_variables(false);

        displacement = displacement_old;

        this->solve();
    }
    catch (...)
    {
        throw;
    }
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::compute_internal_force()
{
    f_int.setZero();

    for (auto const& submesh : fem_mesh.meshes())
    {
        for (std::int64_t element = 0; element < submesh.elements(); ++element)
        {
            auto const [dofs, fe_int] = submesh.internal_force(element);

            f_int(dofs) += fe_int;
        }
    }
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::compute_external_force()
{
    auto const start = std::chrono::high_resolution_clock::now();

    auto const step_time = adaptive_load.step_time();

    f_ext.setZero();

    for (auto const& [name, boundaries] : fem_mesh.nonfollower_boundaries())
    {
        for (auto const& boundary : boundaries.natural_interface())
        {
            std::visit(
                [&](auto const& boundary_mesh) {
                    for (std::int64_t element{0}; element < boundary_mesh.elements(); ++element)
                    {
                        auto const [dofs, fe_ext] = boundary_mesh.external_force(element, step_time);

                        f_ext(dofs) += fe_ext;
                    }
                },
                boundary);
        }
        for (auto const& boundary : boundaries.nodal_interface())
        {
            for (auto dof_index : boundary.dof_view())
            {
                f_ext(dof_index) += boundary.value_view();
            }
        }
    }
    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "External forces assembly took " << elapsed_seconds.count()
              << "s\n";
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::assemble_stiffness()
{
    if (!is_sparsity_computed)
    {
        fem::compute_sparsity_pattern(Kt, fem_mesh);
        is_sparsity_computed = true;
    }

    auto const start = std::chrono::high_resolution_clock::now();

    Kt.coeffs() = 0.0;

    for (auto const& submesh : fem_mesh.meshes())
    {
        tbb::parallel_for(std::int64_t{0}, submesh.elements(), [&](auto const element) {
            auto const [dofs, ke] = submesh.tangent_stiffness(element);

            for (std::int64_t b{0}; b < dofs.size(); b++)
            {
                for (std::int64_t a{0}; a < dofs.size(); a++)
                {
                    Kt.coefficient_update(dofs(a), dofs(b), ke(a, b));
                }
            }
        });
    }

    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Tangent stiffness assembly took "
              << elapsed_seconds.count() << "s\n";
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::enforce_dirichlet_conditions(sparse_matrix& A, vector& b) const
{
    for (auto const& [name, boundaries] : fem_mesh.dirichlet_boundaries())
    {
        for (auto const& dirichlet_boundary : boundaries)
        {
            for (auto const& fixed_dof : dirichlet_boundary.dof_view())
            {
                auto const diagonal_entry = A.coeff(fixed_dof, fixed_dof);

                b(fixed_dof) = 0.0;

                std::vector<std::int32_t> non_zero_visitor;

                // Zero the rows and columns
                for (sparse_matrix::InnerIterator it(A, fixed_dof); it; ++it)
                {
                    // Set the value of the col or row resp. to zero
                    it.valueRef() = 0.0;
                    non_zero_visitor.push_back(A.IsRowMajor ? it.col() : it.row());
                }

                // Zero the row or col respectively
                for (auto const& non_zero : non_zero_visitor)
                {
                    const auto row = A.IsRowMajor ? non_zero : fixed_dof;
                    const auto col = A.IsRowMajor ? fixed_dof : non_zero;

                    A.coeffRef(row, col) = 0.0;
                }
                // Reset the diagonal to the same value to preserve conditioning
                A.coeffRef(fixed_dof, fixed_dof) = diagonal_entry;
            }
        }
    }
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::apply_displacement_boundaries()
{
    Eigen::SparseVector<double> prescribed_increment(displacement.size());

    for (auto const& [name, boundaries] : fem_mesh.dirichlet_boundaries())
    {
        for (auto const& boundary : boundaries)
        {
            auto const delta_u = boundary.value_view(adaptive_load.step_time())
                                 - boundary.value_view(adaptive_load.last_step_time());

            for (auto const& dof : boundary.dof_view())
            {
                prescribed_increment.coeffRef(dof) = delta_u;
            }
        }
    }

    // A sparse matrix - sparse vector multiplication is more efficient for a
    // relatively small vector size with the exception of allocation
    minus_residual -= Kt * prescribed_increment;

    displacement += prescribed_increment;
}

template <class fem_mesh_type>
bool fem_static_matrix<fem_mesh_type>::is_iteration_converged() const
{
    return relative_displacement_norm <= displacement_tolerance
           && relative_force_norm <= residual_tolerance;
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::print_convergence_progress() const
{
    std::cout << std::string(6, ' ') << termcolor::bold;
    if (relative_displacement_norm <= displacement_tolerance)
    {
        std::cout << termcolor::green;
    }
    else
    {
        std::cout << termcolor::yellow;
    }
    std::cout << "Incremental displacement norm " << relative_displacement_norm << "\n"
              << termcolor::reset << std::string(6, ' ');

    if (relative_force_norm <= residual_tolerance)
    {
        std::cout << termcolor::green;
    }
    else
    {
        std::cout << termcolor::yellow;
    }
    std::cout << termcolor::bold << "Residual force norm " << relative_force_norm
              << termcolor::reset << "\n";
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::update_relative_norms()
{
    relative_displacement_norm = delta_d.norm() / displacement.norm();

    relative_force_norm = is_approx(std::max(f_ext.norm(), f_int.norm()), 0.0)
                              ? 1.0
                              : minus_residual.norm() / std::max(f_ext.norm(), f_int.norm());
}

template <class fem_mesh_type>
void fem_static_matrix<fem_mesh_type>::perform_equilibrium_iterations()
{
    displacement = displacement_old;

    // Full Newton-Raphson iteration to solve nonlinear equations
    auto constexpr max_iterations{30};
    auto current_iteration{0};

    while (current_iteration < max_iterations)
    {
        auto const start = std::chrono::high_resolution_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold
                  << "Newton-Raphson iteration " << current_iteration << termcolor::reset
                  << std::endl;

        assemble_stiffness();

        compute_internal_force();

        minus_residual = f_ext - f_int;

        if (current_iteration == 0) apply_displacement_boundaries();

        enforce_dirichlet_conditions(Kt, minus_residual);

        solver->solve(Kt, delta_d, minus_residual);

        displacement += delta_d;

        fem_mesh.update_internal_variables(displacement,
                                           current_iteration == 0 ? adaptive_load.increment() : 0.0);

        update_relative_norms();

        print_convergence_progress();

        auto const end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Equilibrium iteration required "
                  << elapsed_seconds.count() << "s\n";

        if (is_iteration_converged()) break;

        current_iteration++;
    }
    if (current_iteration == max_iterations)
    {
        throw computational_error("Reached Newton-Raphson iteration limit");
    }

    if (current_iteration != max_iterations)
    {
        displacement_old = displacement;

        adaptive_load.update_convergence_state(current_iteration != max_iterations);
        fem_mesh.save_internal_variables(current_iteration != max_iterations);

        fem_mesh.update_internal_forces(f_int);

        io.write(adaptive_load.step(), adaptive_load.time());
    }
}
}
