
#pragma once

/// @file

#include "solver/linear/linear_solver.hpp"
#include "assembler/homogeneous_dirichlet.hpp"
#include "assembler/sparsity_pattern.hpp"
#include "io/file_output.hpp"
#include "io/json.hpp"
#include "solver/adaptive_time_step.hpp"
#include "solver/linear/linear_solver_factory.hpp"

#include <tbb/parallel_for.h>

#include <chrono>
#include <iostream>
#include <variant>

namespace neon
{
/// fem_linear_static_matrix is responsible for assembling the stiffness matrix,
/// solving the linear system of equations and dispatching the postprocessor
/// for a mesh with linear (or infinitesimal) behaviour.  This applies to
/// any mesh that provides the following methods:
/// \p tangent_stiffness
/// \p dirichlet_boundaries
/// and these are used to form the linear system \p K*d=f_ext and this is then
/// passed through the linear solver and the postprocessor.
template <typename fem_mesh_type>
class linear_static_matrix
{
public:
    using mesh_type = fem_mesh_type;

public:
    explicit linear_static_matrix(mesh_type& mesh, json const& simulation_data);

    ~linear_static_matrix() = default;

    void solve();

protected:
    /// Assemble the external contribution vector
    void compute_external_force(double const load_factor = 1.0);

    /// Assemble the tangent stiffness matrix
    void assemble_stiffness();

protected:
    mesh_type& fem_mesh;

    /// Adaptive load stepping
    adaptive_time_step adaptive_load;

    /// Linear solver used to solve the system of equations
    std::unique_ptr<linear_solver> solver;

    /// Stiffness matrix
    sparse_matrix Kt;
    /// Right hand side load vector
    vector f_ext;
    /// Solution vector
    vector d;
    /// Cache the sparsity pattern
    bool is_sparsity_computed{false};
};

template <typename fem_mesh_type>
linear_static_matrix<fem_mesh_type>::linear_static_matrix(mesh_type& fem_mesh,
                                                          json const& simulation_data)
    : fem_mesh(fem_mesh),
      adaptive_load(simulation_data["Time"], fem_mesh.time_history()),
      solver(make_linear_solver(simulation_data["LinearSolver"], fem_mesh.is_symmetric()))
{
    f_ext = d = vector::Zero(fem_mesh.active_dofs());
}

template <typename fem_mesh_type>
void linear_static_matrix<fem_mesh_type>::solve()
{
    fem_mesh.update_internal_variables(d);

    compute_external_force();

    assemble_stiffness();

    apply_dirichlet_conditions(Kt, d, f_ext, fem_mesh);

    solver->solve(Kt, d, f_ext);

    fem_mesh.update_internal_variables(d);

    fem_mesh.write(adaptive_load.step(), adaptive_load.time());
}

template <typename fem_mesh_type>
void linear_static_matrix<fem_mesh_type>::compute_external_force(double)
{
    auto const start = std::chrono::steady_clock::now();

    // auto const step_time = adaptive_load.step_time();

    f_ext.setZero();

    for (auto const& [name, boundaries] : fem_mesh.nonfollower_boundaries())
    {
        // for (auto const& boundary : boundaries.natural_interface())
        // {
        // std::visit(
        //     [&](auto const& boundary_mesh) {
        //         for (std::int64_t element{0}; element < boundary_mesh.elements(); ++element)
        //         {
        //             auto const [dofs, fe_ext] = boundary_mesh.external_force(element, step_time);
        //
        //             f_ext(dofs) += fe_ext;
        //         }
        //     },
        //     boundary);
        // }
        for (auto const& boundary : boundaries.nodal_interface())
        {
            for (auto dof_index : boundary.dof_view())
            {
                f_ext(dof_index) += boundary.value_view();
            }
        }
    }
    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "External forces assembly took " << elapsed_seconds.count()
              << "s\n";
}

template <typename fem_mesh_type>
void linear_static_matrix<fem_mesh_type>::assemble_stiffness()
{
    auto const start = std::chrono::steady_clock::now();

    if (!is_sparsity_computed)
    {
        compute_sparsity_pattern(Kt, fem_mesh);
        is_sparsity_computed = true;
    }

    Kt.coeffs() = 0.0;

    for (auto const& submesh : fem_mesh.meshes())
    {
        tbb::parallel_for(std::int64_t{0}, submesh.elements(), [&](auto const element) {
            auto const [dofs, ke] = submesh.tangent_stiffness(element);

            for (std::int64_t b{0}; b < dofs.size(); b++)
            {
                for (std::int64_t a{0}; a < dofs.size(); a++)
                {
                    Kt.add_to(dofs(a), dofs(b), ke(a, b));
                }
            }
        });
    }

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Tangent stiffness assembly took "
              << elapsed_seconds.count() << "s\n";
}
}
