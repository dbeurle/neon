
#include "femStaticMatrix.hpp"

#include "Exceptions.hpp"
#include "solver/linear/LinearSolverFactory.hpp"

#include <chrono>

#include <json/value.h>
#include <omp.h>
#include <termcolor/termcolor.hpp>

namespace neon::solid
{
femStaticMatrix::femStaticMatrix(femMesh& fem_mesh, Json::Value const& simulation)
    : fem_mesh(fem_mesh),
      file_io(simulation["Name"].asString(), simulation["Visualisation"], fem_mesh),
      adaptive_load(simulation["Time"], fem_mesh.time_history()),
      fint(Vector::Zero(fem_mesh.active_dofs())),
      fext(Vector::Zero(fem_mesh.active_dofs())),
      d(Vector::Zero(fem_mesh.active_dofs())),
      linear_solver(make_linear_solver(simulation["LinearSolver"], fem_mesh.is_symmetric()))
{
    if (!simulation["NonlinearOptions"].isMember("DisplacementTolerance"))
    {
        throw std::runtime_error("DisplacementTolerance not specified in "
                                 "NonlinearOptions");
    }
    if (!simulation["NonlinearOptions"].isMember("ResidualTolerance"))
    {
        throw std::runtime_error("ResidualTolerance not specified in "
                                 "NonlinearOptions");
    }
    residual_tolerance = simulation["NonlinearOptions"]["ResidualTolerance"].asDouble();
    displacement_tolerance = simulation["NonlinearOptions"]["DisplacementTolerance"].asDouble();

    // Perform Newton-Raphson iterations
    std::cout << "\n"
              << std::string(4, ' ') << "Non-linear equation system has " << fem_mesh.active_dofs()
              << " degrees of freedom\n";
}

femStaticMatrix::~femStaticMatrix() = default;

void femStaticMatrix::internal_restart(Json::Value const& solver_data,
                                       Json::Value const& new_increment_data)
{
    adaptive_load.reset(new_increment_data);
    linear_solver = make_linear_solver(solver_data);
}

void femStaticMatrix::compute_sparsity_pattern()
{
    std::vector<Doublet<int>> doublets;
    doublets.reserve(fem_mesh.active_dofs());

    Kt.resize(fem_mesh.active_dofs(), fem_mesh.active_dofs());

    for (auto const& submesh : fem_mesh.meshes())
    {
        // Loop over the elements and add in the non-zero components
        for (auto element = 0; element < submesh.elements(); element++)
        {
            for (auto const& p : submesh.local_dof_list(element))
            {
                for (auto const& q : submesh.local_dof_list(element))
                {
                    doublets.emplace_back(p, q);
                }
            }
        }
    }
    Kt.setFromTriplets(doublets.begin(), doublets.end());
    Kt.finalize();

    is_sparsity_computed = true;
}

void femStaticMatrix::compute_internal_force()
{
    fint.setZero();

    for (auto const& submesh : fem_mesh.meshes())
    {
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            auto const & [ dofs, fe_int ] = submesh.internal_force(element);

            for (auto a = 0; a < fe_int.size(); ++a)
            {
                fint(dofs[a]) += fe_int(a);
            }
        }
    }
}

void femStaticMatrix::compute_external_force(double const step_time)
{
    auto const start = std::chrono::high_resolution_clock::now();

    fext.setZero();

    for (auto const & [ name, nf_loads ] : fem_mesh.nonfollower_load_boundaries())
    {
        for (auto const & [ is_dof_active, boundary_conditions ] : nf_loads.interface())
        {
            if (!is_dof_active) continue;

            for (auto const& boundary_condition : boundary_conditions)
            {
                // clang-format off
                std::visit([&](auto const& mesh) {
                    for (auto element = 0; element < mesh.elements(); ++element)
                    {
                       auto const & [ dofs, fe_ext ] = mesh.external_force(element, step_time);

                       for (auto a = 0; a < fe_ext.size(); ++a)
                       {
                           fext(dofs[a]) += fe_ext(a);
                       }
                    }
                },
                boundary_condition);
                // clang-format on
            }
        }
    }
    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "External forces assembly took " << elapsed_seconds.count()
              << "s\n";
}

void femStaticMatrix::solve()
{
    try
    {
        while (!adaptive_load.is_fully_applied())
        {
            std::cout << "\n"
                      << std::string(4, ' ') << termcolor::magenta << termcolor::bold
                      << "Performing equilibrium iterations for time " << adaptive_load.step_time()
                      << termcolor::reset << std::endl;

            apply_displacement_boundaries();

            compute_external_force(adaptive_load.step_time());

            fem_mesh.update_internal_variables(d, adaptive_load.increment());

            perform_equilibrium_iterations();
        }
    }
    catch (computational_error& comp_error)
    {
        std::cout << std::endl
                  << std::string(6, ' ') << termcolor::bold << termcolor::yellow
                  << comp_error.what() << termcolor::reset << std::endl;

        adaptive_load.update_convergence_state(false);
        fem_mesh.save_internal_variables(false);

        displacement_norm_old = residual_norm_old = 1.0;

        this->solve();
    }
    catch (...)
    {
        throw;
    }
}

void femStaticMatrix::assemble_stiffness()
{
    if (!is_sparsity_computed) compute_sparsity_pattern();

    auto start = std::chrono::high_resolution_clock::now();

    Kt.coeffs() = 0.0;

    for (auto const& submesh : fem_mesh.meshes())
    {
#pragma omp parallel for
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            // auto const[dofs, ke] = submesh.tangent_stiffness(element);
            auto const& tpl = submesh.tangent_stiffness(element);
            auto const& dofs = std::get<0>(tpl);
            auto const& ke = std::get<1>(tpl);

            for (auto b = 0; b < dofs.size(); b++)
            {
                for (auto a = 0; a < dofs.size(); a++)
                {
#pragma omp atomic
                    Kt.coeffRef(dofs[a], dofs[b]) += ke(a, b);
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Tangent stiffness assembly took "
              << elapsed_seconds.count() << "s\n";
}

void femStaticMatrix::enforce_dirichlet_conditions(SparseMatrix& A, Vector& x, Vector& b)
{
    for (auto const & [ name, dirichlet_boundaries ] : fem_mesh.displacement_boundaries())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const& fixed_dof : dirichlet_boundary.dof_view())
            {
                auto const diagonal_entry = A.coeffRef(fixed_dof, fixed_dof);

                x(fixed_dof) = b(fixed_dof) = 0.0;

                std::vector<int> non_zero_visitor;

                // Zero the rows and columns
                for (SparseMatrix::InnerIterator it(A, fixed_dof); it; ++it)
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

void femStaticMatrix::apply_displacement_boundaries()
{
    for (auto const & [ name, dirichlet_boundaries ] : fem_mesh.displacement_boundaries())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const& dof : dirichlet_boundary.dof_view())
            {
                d(dof) = dirichlet_boundary.value_view(adaptive_load.step_time());
            }
        }
    }
}

void femStaticMatrix::perform_equilibrium_iterations()
{
    Vector delta_d = Vector::Zero(fem_mesh.active_dofs());
    Vector d_new = d;

    displacement_norm_old = residual_norm_old = 1.0;

    // Full Newton-Raphson iteration to solve nonlinear equations
    auto constexpr max_iterations{10};
    auto current_iteration{0};

    while (current_iteration < max_iterations)
    {
        auto const start = std::chrono::high_resolution_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold
                  << "Newton-Raphson iteration " << current_iteration << termcolor::reset
                  << std::endl;

        compute_internal_force();

        Vector residual = fint - fext;

        assemble_stiffness();

        enforce_dirichlet_conditions(Kt, delta_d, residual);

        linear_solver->solve(Kt, delta_d, -residual);

        d_new += delta_d;

        fem_mesh.update_internal_variables(d_new);

        update_relative_norms(delta_d.norm(), residual.norm());

        print_convergence_progress();

        auto const end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Equilibrium iteration required "
                  << elapsed_seconds.count() << "s\n";

        if (is_iteration_converged()) break;

        current_iteration++;
    }

    adaptive_load.update_convergence_state(current_iteration != max_iterations);
    fem_mesh.save_internal_variables(current_iteration != max_iterations);

    if (current_iteration != max_iterations)
    {
        d = d_new;
        file_io.write(adaptive_load.step(), adaptive_load.time());
    }
}

void femStaticMatrix::update_relative_norms(double const delta_d_norm, double const residual_norm)
{
    relative_displacement_norm = delta_d_norm / displacement_norm_old;
    relative_force_norm = residual_norm / residual_norm_old;

    displacement_norm_old = delta_d_norm;
    residual_norm_old = residual_norm;
}

void femStaticMatrix::print_convergence_progress() const
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
}
