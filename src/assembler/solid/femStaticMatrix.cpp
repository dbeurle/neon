
#include "femStaticMatrix.hpp"

#include "solver/linear/LinearSolverFactory.hpp"

#include <chrono>

#include <omp.h>
#include <termcolor/termcolor.hpp>

namespace neon::solid
{
femStaticMatrix::femStaticMatrix(femMesh& fem_mesh,
                                 Json::Value const& solver_data,
                                 Json::Value const& increment_data)
    : fem_mesh(fem_mesh),
      adaptive_load(increment_data),
      fint(Vector::Zero(fem_mesh.active_dofs())),
      d(Vector::Zero(fem_mesh.active_dofs())),
      linear_solver(make_linear_solver(solver_data))
{
}

femStaticMatrix::~femStaticMatrix() = default;

void femStaticMatrix::continuation(Json::Value const& new_increment_data)
{
    adaptive_load.reset(new_increment_data);
}

void femStaticMatrix::compute_sparsity_pattern()
{
    // auto start = std::chrono::high_resolution_clock::now();

    std::vector<Doublet> doublets;
    doublets.reserve(fem_mesh.active_dofs());

    Kt.resize(fem_mesh.active_dofs(), fem_mesh.active_dofs());

    for (const auto& submesh : fem_mesh.meshes())
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

    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end - start;

    // std::cout << "  Sparsity pattern with " << Kt.nonZeros() << " non-zeros took "
    //           << elapsed_seconds.count() << "s\n";
}

void femStaticMatrix::solve()
{
    // Perform Newton-Raphson iterations
    std::cout << std::string(4, ' ') << "Non-linear equation system has " << fem_mesh.active_dofs()
              << " degrees of freedom\n";

    if (adaptive_load.step() == 0) fem_mesh.write(0, 0.0);

    while (!adaptive_load.is_fully_applied())
    {
        std::cout << termcolor::magenta << termcolor::bold << "\n"
                  << std::string(4, ' ') << "Performing equilibrium iterations for time "
                  << adaptive_load.step_time() << termcolor::reset << std::endl;

        apply_displacement_boundaries();

        fem_mesh.update_internal_variables(d, adaptive_load.increment());

        perform_equilibrium_iterations();
    }
}

void femStaticMatrix::assemble_stiffness()
{
    using ranges::accumulate;

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

void femStaticMatrix::compute_internal_force()
{
    // auto start = std::chrono::high_resolution_clock::now();

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

    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end - start;

    // std::cout << "  Assembly of internal forces took " << elapsed_seconds.count() << "s\n";
}

void femStaticMatrix::enforce_dirichlet_conditions(SparseMatrix& A, Vector& x, Vector& b)
{
    // auto start = std::chrono::high_resolution_clock::now();

    for (auto const & [ name, dirichlet_boundaries ] : fem_mesh.dirichlet_boundary_map())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const& fixed_dof : dirichlet_boundary.dof_view())
            {
                auto const diagonal_entry = Kt.coeffRef(fixed_dof, fixed_dof);

                x(fixed_dof) = b(fixed_dof) = 0.0;

                std::vector<int> non_zero_visitor;

                // Zero the rows and columns
                for (SparseMatrix::InnerIterator it(Kt, fixed_dof); it; ++it)
                {
                    // Set the value of the col or row resp. to zero
                    it.valueRef() = 0.0;
                    non_zero_visitor.push_back(Kt.IsRowMajor ? it.col() : it.row());
                }

                // Zero the row or col respectively
                for (auto const& non_zero : non_zero_visitor)
                {
                    const auto row = Kt.IsRowMajor ? non_zero : fixed_dof;
                    const auto col = Kt.IsRowMajor ? fixed_dof : non_zero;

                    Kt.coeffRef(row, col) = 0.0;
                }
                // Reset the diagonal to the same value to preserve conditioning
                Kt.coeffRef(fixed_dof, fixed_dof) = diagonal_entry;
            }
        }
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end - start;

    // std::cout << "  Dirichlet conditions enforced in " << elapsed_seconds.count() << "s\n";
}

void femStaticMatrix::apply_displacement_boundaries()
{
    // auto start = std::chrono::high_resolution_clock::now();

    for (auto const & [ name, dirichlet_boundaries ] : fem_mesh.dirichlet_boundary_map())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const& dof : dirichlet_boundary.dof_view())
            {
                d(dof) = dirichlet_boundary.value_view(adaptive_load.load_factor());
            }
        }
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end - start;

    // std::cout << "  Displacements applied in " << elapsed_seconds.count() << "s\n";
}

void femStaticMatrix::perform_equilibrium_iterations()
{
    Vector delta_d = Vector::Zero(fem_mesh.active_dofs());

    // Full Newton-Raphson iteration to solve nonlinear equations
    auto constexpr max_iterations = 10;
    auto current_iteration = 0;

    while (current_iteration < max_iterations)
    {
        auto start = std::chrono::high_resolution_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold
                  << "Newton-Raphson iteration " << current_iteration << termcolor::reset
                  << std::endl;

        compute_internal_force();

        Vector residual = fint; // - fext

        assemble_stiffness();

        enforce_dirichlet_conditions(Kt, delta_d, residual);

        linear_solver->solve(Kt, delta_d, -residual);

        d += delta_d;

        fem_mesh.update_internal_variables(d, adaptive_load.increment());

        print_convergence_progress(delta_d.norm(), residual.norm());

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Equilibrium iteration required "
                  << elapsed_seconds.count() << "s\n";

        if (is_converged(delta_d.norm(), residual.norm()))
        {
            break;
        }
        current_iteration++;
    }
    adaptive_load.update_convergence_state(current_iteration != max_iterations);
    fem_mesh.save_internal_variables(current_iteration != max_iterations);

    std::cout << "\n"
              << std::string(4, ' ') << "Writing solution to file for step " << adaptive_load.step()
              << "\n";

    if (current_iteration != max_iterations)
        fem_mesh.write(adaptive_load.step(), adaptive_load.time());
}

void femStaticMatrix::print_convergence_progress(double const delta_d_norm,
                                                 double const residual_norm) const
{
    std::cout << termcolor::bold;
    if (delta_d_norm <= displacement_tolerance)
    {
        std::cout << termcolor::green;
    }
    else
    {
        std::cout << termcolor::yellow;
    }
    std::cout << std::string(6, ' ') << "Incremental displacement norm " << delta_d_norm << "\n";
    std::cout << termcolor::reset << termcolor::bold;

    if (residual_norm <= residual_tolerance)
    {
        std::cout << termcolor::green;
    }
    else
    {
        std::cout << termcolor::yellow;
    }
    std::cout << std::string(6, ' ') << "Residual force norm " << residual_norm << termcolor::reset
              << "\n";
}
}
