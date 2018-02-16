
#include "fem_static_matrix.hpp"

#include "Exceptions.hpp"
#include "solver/linear/LinearSolverFactory.hpp"

#include "io/json.hpp"
#include <chrono>
#include <omp.h>

namespace neon::diffusion
{
fem_static_matrix::fem_static_matrix(fem_mesh& mesh, json const& simulation_data)
    : mesh(mesh),
      f(vector::Zero(mesh.active_dofs())),
      d(vector::Zero(mesh.active_dofs())),
      file_io(simulation_data["Name"].get<std::string>(), simulation_data["Visualisation"], mesh),
      linear_solver(make_linear_solver(simulation_data["LinearSolver"]))
{
}

fem_static_matrix::~fem_static_matrix() = default;

void fem_static_matrix::compute_sparsity_pattern()
{
    std::vector<Doublet<int>> doublets;
    doublets.reserve(mesh.active_dofs());

    K.resize(mesh.active_dofs(), mesh.active_dofs());

    for (auto const& submesh : mesh.meshes())
    {
        // Loop over the elements and add in the non-zero components
        for (std::size_t element{0}; element < submesh.elements(); element++)
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
    K.setFromTriplets(doublets.begin(), doublets.end());
    is_sparsity_computed = true;
}

void fem_static_matrix::compute_external_force(double const load_factor)
{
    auto const start = std::chrono::high_resolution_clock::now();

    f.setZero();

    for (auto const& [name, surfaces] : mesh.surface_boundaries())
    {
        for (auto const& surface : surfaces)
        {
            for (auto const& mesh : surface.load_interface())
            {
                for (std::size_t element{0}; element < mesh.elements(); ++element)
                {
                    auto const& [dofs, fe] = mesh.external_force(element, load_factor);

                    for (auto a = 0l; a < fe.size(); ++a)
                    {
                        f(dofs[a]) += fe(a);
                    }
                }
            }
            for (auto const& mesh : surface.stiffness_load_interface())
            {
                for (std::size_t element{0}; element < mesh.elements(); ++element)
                {
                    auto const& [dofs, fe] = mesh.external_force(element, load_factor);
                    auto const& [_, ke] = mesh.external_stiffness(element, load_factor);

                    for (auto a = 0l; a < fe.size(); ++a)
                    {
                        f(dofs[a]) += fe(a);

                        for (auto b{0l}; b < fe.size(); ++b)
                        {
                            K.coeffRef(dofs[a], dofs[b]) += ke(a, b);
                        }
                    }
                }
            }
        }
    }
    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "External forces assembly took " << elapsed_seconds.count()
              << "s\n";
}

void fem_static_matrix::solve()
{
    std::cout << std::string(4, ' ') << "Linear equation system has " << mesh.active_dofs()
              << " degrees of freedom\n";

    assemble_stiffness();

    compute_external_force();

    apply_dirichlet_conditions(K, d, f);

    linear_solver->solve(K, d, f);

    file_io.write(0, 0.0, d);
}

void fem_static_matrix::assemble_stiffness()
{
    if (!is_sparsity_computed) compute_sparsity_pattern();

    auto const start = std::chrono::high_resolution_clock::now();

    K.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
#pragma omp parallel for
        for (std::size_t element = 0; element < submesh.elements(); ++element)
        {
            // auto const[dofs, ke] = submesh.tangent_stiffness(element);
            auto const& tpl = submesh.tangent_stiffness(element);
            auto const& dofs = std::get<0>(tpl);
            auto const& ke = std::get<1>(tpl);

            for (std::size_t b = 0; b < dofs.size(); b++)
            {
                for (std::size_t a = 0; a < dofs.size(); a++)
                {
#pragma omp atomic
                    K.coeffRef(dofs[a], dofs[b]) += ke(a, b);
                }
            }
        }
    }

    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Stiffness assembly took " << elapsed_seconds.count()
              << "s\n";
}

void fem_static_matrix::apply_dirichlet_conditions(sparse_matrix& A, vector& x, vector& b)
{
    for (auto const& [name, dirichlet_boundaries] : mesh.dirichlet_boundaries())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const& fixed_dof : dirichlet_boundary.dof_view())
            {
                auto const diagonal_entry = A.coeffRef(fixed_dof, fixed_dof);

                x(fixed_dof) = dirichlet_boundary.value_view();

                std::vector<int> non_zero_visitor;

                for (sparse_matrix::InnerIterator it(A, fixed_dof); it; ++it)
                {
                    if (!A.IsRowMajor) b(it.row()) -= it.valueRef() * x(fixed_dof);

                    it.valueRef() = 0.0;

                    non_zero_visitor.push_back(A.IsRowMajor ? it.col() : it.row());
                }

                for (auto const& non_zero : non_zero_visitor)
                {
                    auto const row = A.IsRowMajor ? non_zero : fixed_dof;
                    auto const col = A.IsRowMajor ? fixed_dof : non_zero;

                    if (A.IsRowMajor) b(row) -= A.coeffRef(row, col) * x(fixed_dof);

                    A.coeffRef(row, col) = 0.0;
                }

                // Reset the diagonal to the same value to preserve condition number
                A.coeffRef(fixed_dof, fixed_dof) = diagonal_entry;
                b(fixed_dof) = diagonal_entry * x(fixed_dof);
            }
        }
    }
}
}
