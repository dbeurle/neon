
#pragma once

/// @file

#include "static_matrix.hpp"

namespace neon::mechanics
{
/// Generic latin matrix designed for solid mechanics problems using the
/// LATIN method for the solution of nonlinear equations.
/// This class is responsible for the assembly of the process stiffness matrix,
/// residual force vector and solution of the incremental displacement.
template <class MeshType>
class latin_matrix : public static_matrix<MeshType>
{
public:
    using mesh_type = MeshType;
    using base_type = static_matrix<mesh_type>;

private:
    // use base class variables and methods
    using base_type::adaptive_load;
    using base_type::apply_displacement_boundaries;
    using base_type::assemble_stiffness;
    using base_type::compute_external_force;
    using base_type::compute_internal_force;
    using base_type::delta_d;
    using base_type::displacement;
    using base_type::displacement_old;
    using base_type::f_ext;
    using base_type::f_int;
    using base_type::is_iteration_converged;
    using base_type::Kt;
    using base_type::maximum_iterations;
    using base_type::mesh;
    using base_type::minus_residual;
    using base_type::norm_initial_residual;
    using base_type::print_convergence_progress;
    using base_type::solver;
    using base_type::update_relative_norms;

private:
    /// LATIN residual vector
    vector latin_residual;

    /// LATIN search direction scaling factor \f$ \alpha \ \mathbb{C} \f$
    double latin_search_direction = 0.0;

public:
    explicit latin_matrix(mesh_type& mesh, json const& simulation);

    /// Solve the nonlinear system of equations
    void solve();

private:
    void perform_equilibrium_iterations();

    /// Gathers latin internal force vector using current and old Cauchy stresses
    void compute_incremental_latin_internal_force();

    /// Apply dirichlet conditions to the system defined by A, x, and (b or c).
    /// This method sets the incremental displacements to zero for the given
    /// load increment such that incremental displacements are zero
    void enforce_dirichlet_conditions(sparse_matrix& A, vector& b, vector& c) const;
};

template <class MeshType>
latin_matrix<MeshType>::latin_matrix(mesh_type& mesh, json const& simulation)
    : base_type::static_matrix(mesh, simulation)
{
    auto const& nonlinear_options = simulation["nonlinear_options"];

    if (nonlinear_options.find("latin_search_direction") == nonlinear_options.end())
    {
        throw std::domain_error("latin_search_direction not specified in nonlinear_options");
    }

    latin_search_direction = nonlinear_options["latin_search_direction"];

    latin_residual = vector::Zero(mesh.active_dofs());
}

template <class MeshType>
void latin_matrix<MeshType>::compute_incremental_latin_internal_force()
{
    latin_residual.setZero();

    for (auto const& submesh : mesh.meshes())
    {
        for (std::int64_t element{0}; element < submesh.elements(); ++element)
        {
            auto const
                & [dofs, fe_int] = submesh.incremental_latin_internal_force(element,
                                                                            latin_search_direction);

            latin_residual(dofs) += fe_int;
        }
    }
}

template <class MeshType>
void latin_matrix<MeshType>::enforce_dirichlet_conditions(sparse_matrix& A, vector& b, vector& c) const
{
    for (auto const& [name, boundaries] : mesh.dirichlet_boundaries())
    {
        for (auto const& boundary : boundaries)
        {
            if (boundary.is_not_active(adaptive_load.step_time()))
            {
                continue;
            }

            for (auto const& fixed_dof : boundary.dof_view())
            {
                auto const diagonal_entry = A.coeff(fixed_dof, fixed_dof);

                b(fixed_dof) = 0.0;

                c(fixed_dof) = 0.0;

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

template <class MeshType>
void latin_matrix<MeshType>::solve()
{
    try
    {
        // Initialise the mesh with zero displacements
        mesh.update_internal_variables(displacement);
        mesh.update_internal_forces(f_int);

        mesh.write(adaptive_load.step(), adaptive_load.time());

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
        mesh.save_internal_variables(false);

        displacement = displacement_old;

        this->solve();
    }
    catch (...)
    {
        throw;
    }
}

template <class MeshType>
void latin_matrix<MeshType>::perform_equilibrium_iterations()
{
    displacement = displacement_old;

    mesh.update_internal_variables(displacement, adaptive_load.increment());

    // Full LATIN iteration to solve nonlinear equations
    auto current_iteration{0};
    while (current_iteration < maximum_iterations)
    {
        auto const start = std::chrono::steady_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold << "LATIN iteration "
                  << current_iteration << termcolor::reset << "\n";

        assemble_stiffness();

        compute_internal_force();

        compute_incremental_latin_internal_force();

        minus_residual = f_ext - f_int;

        if (current_iteration == 0)
        {
            apply_displacement_boundaries();
            norm_initial_residual = minus_residual.norm();

            // start with an elastic initialisation
            latin_residual = minus_residual;
        }

        // TODO: the convergence may be measure by the minus_residual = latin_residual. However, this
        // will not ensure the balance of forces anymore. Also, the stifness matrix may be scaled by
        // `latin_search_direction` but this did not improve the convergence of the incremental LATIN scheme
        enforce_dirichlet_conditions(Kt, minus_residual, latin_residual);

        solver->solve(Kt, delta_d, latin_residual);

        displacement += delta_d;

        mesh.update_internal_variables(displacement, 0.0);

        update_relative_norms();

        print_convergence_progress();

        auto const end = std::chrono::steady_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Equilibrium iteration required "
                  << elapsed_seconds.count() << "s\n";

        if (is_iteration_converged()) break;

        current_iteration++;
    }
    if (current_iteration == maximum_iterations)
    {
        throw computational_error("Reached LATIN iteration limit");
    }

    if (current_iteration != maximum_iterations)
    {
        displacement_old = displacement;

        adaptive_load.update_convergence_state(current_iteration != maximum_iterations);
        mesh.save_internal_variables(current_iteration != maximum_iterations);

        mesh.update_internal_forces(f_int);

        mesh.write(adaptive_load.step(), adaptive_load.time());
    }
}

}
