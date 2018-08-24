
#include "static_matrix.hpp"

#include "exceptions.hpp"
#include "solver/linear/linear_solver_factory.hpp"
#include "assembler/homogeneous_dirichlet.hpp"
#include "numeric/doublet.hpp"
#include "io/json.hpp"

#include <tbb/parallel_for.h>

#include <chrono>

namespace neon::diffusion
{
static_matrix::static_matrix(mesh_type& mesh, json const& simulation_data)
    : mesh(mesh),
      f(vector::Zero(mesh.active_dofs())),
      d(vector::Zero(mesh.active_dofs())),
      solver(make_linear_solver(simulation_data["LinearSolver"]))
{
}

static_matrix::~static_matrix() = default;

void static_matrix::compute_sparsity_pattern()
{
    std::vector<doublet<std::int32_t>> doublets;
    doublets.reserve(mesh.active_dofs());

    K.resize(mesh.active_dofs(), mesh.active_dofs());

    for (auto const& submesh : mesh.meshes())
    {
        // Loop over the elements and add in the non-zero components
        for (std::int64_t element{0}; element < submesh.elements(); element++)
        {
            auto const dof_view = submesh.local_dof_view(element);

            for (std::int64_t p{0}; p < dof_view.size(); ++p)
            {
                for (std::int64_t q{0}; q < dof_view.size(); ++q)
                {
                    doublets.emplace_back(dof_view(p), dof_view(q));
                }
            }
        }
    }
    K.setFromTriplets(begin(doublets), end(doublets));

    is_sparsity_computed = true;
}

void static_matrix::compute_external_force(double const load_factor)
{
    auto const start = std::chrono::steady_clock::now();

    f.setZero();

    for (auto const& [name, surfaces] : mesh.surface_boundaries())
    {
        for (auto const& surface : surfaces)
        {
            for (auto const& mesh : surface.load_interface())
            {
                for (std::int64_t element{0}; element < mesh.elements(); ++element)
                {
                    auto const [dof_view, fe] = mesh.external_force(element, load_factor);

                    f(dof_view) += fe;
                }
            }
            for (auto const& mesh : surface.stiffness_load_interface())
            {
                for (std::int64_t element{0}; element < mesh.elements(); ++element)
                {
                    auto const [dof_view, fe] = mesh.external_force(element, load_factor);
                    auto const [_, ke] = mesh.external_stiffness(element, load_factor);

                    f(dof_view) += fe;

                    for (std::int64_t a{0}; a < fe.size(); ++a)
                    {
                        for (std::int64_t b{0}; b < fe.size(); ++b)
                        {
                            K.coefficient_update(dof_view(a), dof_view(b), ke(a, b));
                        }
                    }
                }
            }
        }
    }
    auto const end = std::chrono::steady_clock::now();

    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "External forces assembly took " << elapsed_seconds.count()
              << "s\n";
}

void static_matrix::solve()
{
    mesh.update_internal_variables(d);

    mesh.write(0, 0.0);

    std::cout << std::string(4, ' ') << "Linear equation system has " << mesh.active_dofs()
              << " degrees of freedom\n";

    assemble_stiffness();

    compute_external_force();

    fem::apply_dirichlet_conditions(K, d, f, mesh);

    solver->solve(K, d, f);

    mesh.update_internal_variables(d);

    mesh.write(1, 1.0);
}

void static_matrix::assemble_stiffness()
{
    if (!is_sparsity_computed) compute_sparsity_pattern();

    auto const start = std::chrono::steady_clock::now();

    K.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
        tbb::parallel_for(std::int64_t{0}, submesh.elements(), [&](auto const element) {
            auto const& [dof_view, local_tangent] = submesh.tangent_stiffness(element);

            for (std::int64_t a{0}; a < dof_view.size(); a++)
            {
                for (std::int64_t b{0}; b < dof_view.size(); b++)
                {
                    K.coefficient_update(dof_view(a), dof_view(b), local_tangent(a, b));
                }
            }
        });
    }

    auto const end = std::chrono::steady_clock::now();

    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Stiffness assembly took " << elapsed_seconds.count()
              << "s\n";
}
}
