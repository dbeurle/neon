
#include "femStaticMatrix.hpp"

#include "Exceptions.hpp"
#include "numeric/float_compare.hpp"
#include "solver/linear/LinearSolverFactory.hpp"

#include <cfenv>
#include <chrono>

#include "io/json.hpp"
#include <omp.h>
#include <termcolor/termcolor.hpp>

namespace neon::mechanical::solid
{
femStaticMatrix::femStaticMatrix(fem_mesh& mesh, json const& simulation)
    : mesh(mesh),
      io(simulation["Name"].get<std::string>(), simulation["Visualisation"], mesh),
      adaptive_load(simulation["Time"], mesh.time_history()),
      fint(vector::Zero(mesh.active_dofs())),
      fext(vector::Zero(mesh.active_dofs())),
      displacement(vector::Zero(mesh.active_dofs())),
      displacement_old(vector::Zero(mesh.active_dofs())),
      delta_d(vector::Zero(mesh.active_dofs())),
      linear_solver(make_linear_solver(simulation["LinearSolver"], mesh.is_symmetric()))
{
    if (!simulation["NonlinearOptions"].count("DisplacementTolerance"))
    {
        throw std::runtime_error("DisplacementTolerance not specified in "
                                 "NonlinearOptions");
    }
    if (!simulation["NonlinearOptions"].count("ResidualTolerance"))
    {
        throw std::runtime_error("ResidualTolerance not specified in "
                                 "NonlinearOptions");
    }
    residual_tolerance = simulation["NonlinearOptions"]["ResidualTolerance"];
    displacement_tolerance = simulation["NonlinearOptions"]["DisplacementTolerance"];

    // Perform Newton-Raphson iterations
    std::cout << "\n"
              << std::string(4, ' ') << "Non-linear equation system has " << mesh.active_dofs()
              << " degrees of freedom\n";
}

femStaticMatrix::~femStaticMatrix() = default;

void femStaticMatrix::internal_restart(json const& solver_data, json const& new_increment_data)
{
    adaptive_load.reset(new_increment_data);
    linear_solver = make_linear_solver(solver_data);
}

void femStaticMatrix::compute_sparsity_pattern()
{
    std::vector<Doublet<std::int32_t>> doublets;
    doublets.reserve(mesh.active_dofs());

    Kt.resize(mesh.active_dofs(), mesh.active_dofs());

    for (auto const& submesh : mesh.meshes())
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

    std::feclearexcept(FE_ALL_EXCEPT);

    for (auto const& submesh : mesh.meshes())
    {
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            auto const& [dofs, fe_int] = submesh.internal_force(element);

            for (auto a = 0; a < fe_int.size(); ++a)
            {
                fint(dofs[a]) += fe_int(a);
            }
        }
    }
    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }
}

void femStaticMatrix::compute_external_force()
{
    auto const start = std::chrono::high_resolution_clock::now();

    auto const step_time = adaptive_load.step_time();

    std::feclearexcept(FE_ALL_EXCEPT);

    fext.setZero();

    for (auto const& [name, nf_loads] : mesh.nonfollower_load_boundaries())
    {
        for (auto const& [is_dof_active, boundary_conditions] : nf_loads.interface())
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
    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
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
        // Initialise the mesh with zero displacements
        mesh.update_internal_variables(displacement);

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

void femStaticMatrix::assemble_stiffness()
{
    if (!is_sparsity_computed) compute_sparsity_pattern();

    auto const start = std::chrono::high_resolution_clock::now();

    Kt.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
#pragma omp parallel for
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            // auto const[dofs, ke] = submesh.tangent_stiffness(element);
            auto const& tpl = submesh.tangent_stiffness(element);
            auto const dofs = std::get<0>(tpl);
            auto const ke = std::get<1>(tpl);

            for (auto a = 0; a < dofs.size(); a++)
            {
                for (auto b = 0; b < dofs.size(); b++)
                {
#pragma omp atomic
                    Kt.coeffRef(dofs[a], dofs[b]) += ke(a, b);
                }
            }
        }
    }

    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Tangent stiffness assembly took "
              << elapsed_seconds.count() << "s\n";
}

void femStaticMatrix::enforce_dirichlet_conditions(sparse_matrix& A, vector& b) const
{
    for (auto const& [name, boundaries] : mesh.displacement_boundaries())
    {
        for (auto const& dirichlet_boundary : boundaries)
        {
            for (auto const& fixed_dof : dirichlet_boundary.dof_view())
            {
                auto const diagonal_entry = A.coeff(fixed_dof, fixed_dof);

                b(fixed_dof) = 0.0;

                std::vector<int> non_zero_visitor;

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

void femStaticMatrix::apply_displacement_boundaries()
{
    Eigen::SparseVector<double> prescribed_increment(displacement.size());

    for (auto const& [name, boundaries] : mesh.displacement_boundaries())
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

void femStaticMatrix::perform_equilibrium_iterations()
{
    displacement = displacement_old;

    // Full Newton-Raphson iteration to solve nonlinear equations
    auto constexpr max_iterations{10};
    auto current_iteration{0};

    while (current_iteration < max_iterations)
    {
        auto const start = std::chrono::high_resolution_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold
                  << "Newton-Raphson iteration " << current_iteration << termcolor::reset
                  << std::endl;

        assemble_stiffness();

        compute_internal_force();

        minus_residual = fext - fint;

        if (current_iteration == 0) apply_displacement_boundaries();

        enforce_dirichlet_conditions(Kt, minus_residual);

        linear_solver->solve(Kt, delta_d, minus_residual);

        displacement += delta_d;

        mesh.update_internal_variables(displacement,
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
    if (current_iteration != max_iterations)
    {
        adaptive_load.update_convergence_state(current_iteration != max_iterations);
        mesh.save_internal_variables(current_iteration != max_iterations);

        displacement_old = displacement;
        io.write(adaptive_load.step(), adaptive_load.time());
    }
    else
    {
        throw computational_error("Reached Newton-Raphson iteration limit");
    }
}

void femStaticMatrix::update_relative_norms()
{
    relative_displacement_norm = delta_d.norm() / displacement.norm();

    relative_force_norm = is_approx(std::max(fext.norm(), fint.norm()), 0.0)
                              ? 1.0
                              : minus_residual.norm() / std::max(fext.norm(), fint.norm());
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
