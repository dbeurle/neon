
#include "femMatrix.hpp"

#include "solver/linear/LinearSolverFactory.hpp"

#include <boost/timer/timer.hpp>
#include <iostream>

namespace neon::solid
{
femMatrix::femMatrix(femMesh& fem_mesh, Json::Value const& solver_data)
    : fem_mesh(fem_mesh),
      fint(fem_mesh.active_dofs()),
      linear_solver(LinearSolverFactory::make(solver_data))
{
}

femMatrix::~femMatrix() = default;

void femMatrix::compute_sparsity_pattern()
{
    boost::timer::cpu_timer timer;

    std::vector<Doublet> doublets;

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
                    doublets.push_back({p, q});
                }
            }
        }
    }

    Kt.setFromTriplets(doublets.begin(), doublets.end());
    Kt.finalize();

    is_sparsity_computed = true;

    std::cout << "  Sparsity pattern with " << Kt.nonZeros() << " non-zeros took" << timer.format();
}

void femMatrix::solve()
{
    // Perform Newton-Raphson iterations
    std::cout << "Solving " << fem_mesh.active_dofs() << " non-linear equations\n";

    Vector delta_d = Vector::Zero(fem_mesh.active_dofs());

    // Iteration and tolerance specifications
    auto max_iterations = 15, current_iteration = 0;
    auto tolerance = 1.0e-5;

    apply_displacements();

    // Full Newton-Raphson iteration to solve nonlinear equations
    while (current_iteration < max_iterations)
    {
        std::cout << "----------------------------------\n";
        std::cout << "    Newton-Raphson iteration " << current_iteration << "\n";
        std::cout << "----------------------------------\n";

        compute_internal_force();

        Vector residual = -fint;

        assemble_stiffness();

        enforce_dirichlet_conditions(Kt, delta_d, residual);

        std::cout << "  Displacement norm " << delta_d.norm() << ", residual norm "
                  << residual.norm() << "\n";

        if (delta_d.norm() < tolerance && residual.norm() < tolerance)
        {
            std::cout << "Nonlinear iterations converged!\n";
            break;
        }

        linear_solver->solve(Kt, delta_d, residual);

        fem_mesh.update_internal_variables(delta_d);

        fem_mesh.write();

        current_iteration++;
    }

    if (current_iteration == max_iterations)
    {
        std::cout << "Newton-Raphson iterations failed to converged.  Aborting.\n";
        std::abort();
    }
}

void femMatrix::assemble_stiffness()
{
    if (!is_sparsity_computed) compute_sparsity_pattern();

    boost::timer::cpu_timer timer;

    Kt.coeffs() = 0.0;

    for (const auto& submesh : fem_mesh.meshes())
    {
        // #pragma omp parallel for
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            auto const[dofs, ke] = submesh.tangent_stiffness(element);

            // Gather and add to the global coefficient matrix
            for (auto b = 0; b < dofs.size(); b++)
            {
                for (auto a = 0; a < dofs.size(); a++)
                {
                    // #pragma omp atomic update
                    Kt.coeffRef(dofs[a], dofs[b]) += ke(a, b);
                }
            }
        }
    }

    std::cout << "  Assembly of "
              << ranges::accumulate(fem_mesh.meshes(),
                                    0l,
                                    [](auto sum, auto const& submesh) {
                                        return sum + submesh.elements();
                                    })
              << " elements with " << Kt.nonZeros() << " insertions took" << timer.format();
}

void femMatrix::compute_internal_force()
{
    boost::timer::cpu_timer timer;

    fint.setZero();

    for (const auto& submesh : fem_mesh.meshes())
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
    std::cout << "  Assembly of internal forces took" << timer.format();
}

void femMatrix::enforce_dirichlet_conditions(SparseMatrix& A, Vector& x, Vector& b)
{
    boost::timer::cpu_timer timer;

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
                for (const auto& non_zero : non_zero_visitor)
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
    std::cout << "  Dirichlet conditions enforced in" << timer.format();
}

void femMatrix::apply_displacements(double const load_factor)
{
    boost::timer::cpu_timer timer;

    Vector disp = Vector::Zero(fem_mesh.active_dofs());

    for (const auto & [ name, dirichlet_boundaries ] : fem_mesh.dirichlet_boundary_map())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const& dof : dirichlet_boundary.dof_view())
            {
                disp(dof) = load_factor * dirichlet_boundary.value_view();
            }
        }
    }
    fem_mesh.update_internal_variables(disp);
    std::cout << "  Displacements applied in" << timer.format();
}
}
