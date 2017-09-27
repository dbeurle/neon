
#include "femDynamicMatrix.hpp"

#include "solver/linear/LinearSolver.hpp"

#include <chrono>
#include <json/value.h>
#include <termcolor/termcolor.hpp>

namespace neon::diffusion
{
femDynamicMatrix::femDynamicMatrix(femMesh& fem_mesh,
                                   Json::Value const& simulation_data,
                                   FileIO&& file_io)
    : femStaticMatrix(fem_mesh, simulation_data, std::move(file_io)),
      time_solver(simulation_data["Time"])
{
}

void femDynamicMatrix::solve()
{
    // Perform time dependent solution
    std::cout << "Solving " << fem_mesh.active_dofs() << " degrees of freedom\n";

    assemble_mass();

    std::cout << M << std::endl;

    compute_external_force();

    assemble_stiffness();

    std::cout << std::endl << K << std::endl;

    SparseMatrix A = M - time_solver.current_time_step_size() * K;

    while (time_solver.loop())
    {
        std::cout << "Performing time step\n";

        apply_dirichlet_conditions(A, d, f);

        linear_solver->solve(A, d, f);

        // perform_equilibrium_iterations();

        // fem_mesh.write();
    }
    std::cout << "Solver routine finished\n";
}

void femDynamicMatrix::assemble_mass()
{
    auto start = std::chrono::high_resolution_clock::now();

    M.resize(fem_mesh.active_dofs(), fem_mesh.active_dofs());

    std::vector<Doublet<int>> doublets;
    doublets.reserve(fem_mesh.active_dofs());

    for (auto const& submesh : fem_mesh.meshes())
    {
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
    M.setFromTriplets(doublets.begin(), doublets.end());

    for (auto const& submesh : fem_mesh.meshes())
    {
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            auto const[dofs, m] = submesh.consistent_mass(element);

            for (auto b = 0; b < dofs.size(); b++)
            {
                for (auto a = 0; a < dofs.size(); a++)
                {
                    M.coeffRef(dofs[a], dofs[b]) += m(a, b);
                }
            }
        }
    }

    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(4, ' ') << "Assembly of mass matrix took " << elapsed_seconds.count()
              << "s\n";
}

// void femDynamicMatrix::perform_equilibrium_iterations()
// {
//     Vector delta_d = Vector::Zero(fem_mesh.active_dofs());
//
//     // Full Newton-Raphson iteration to solve nonlinear equations
//     auto constexpr max_iterations = 20;
//     auto current_iteration = 0;
//     auto constexpr tolerance = 1.0e-5;
//
//     while (current_iteration < max_iterations)
//     {
//         std::cout << "\n  Newton-Raphson iteration " << current_iteration << "\n\n";
//
//         compute_internal_force();
//
//         a = newmark.accelerations(a, v, d);
//         v = newmark.velocities(a, v);
//
//         std::cout << "\nAcceleration\n\n" << a << std::endl;
//         std::cout << "\nVelocity\n\n" << v << std::endl;
//         std::cout << "\nDisplacement\n\n" << d << std::endl;
//
//         auto f = -fint;
//
//         Vector residual = M.cwiseProduct(a) - f;
//
//         // Build A = 1 / (β Δt^2) * M + K_int - K_ext
//         assemble_stiffness();
//         Kt.diagonal() += newmark.mass_scaling_factor() * M.diagonal();
//
//         enforce_dirichlet_conditions(Kt, delta_d, residual);
//
//         std::cout << "\nResidual\n\n" << residual << std::endl;
//
//         // if (auto norm = delta_d.norm(); norm > tolerance)
//         //     std::cout <<
//
//         std::cout << "    Displacement norm " << delta_d.norm() << ", residual norm "
//                   << residual.norm() << "\n";
//
//         if (delta_d.norm() < tolerance && residual.norm() < tolerance) break;
//
//         linear_solver->solve(Kt, delta_d, -residual);
//
//         d += delta_d;
//
//         fem_mesh.update_internal_variables(d, newmark.time_step_size());
//
//         // fem_mesh.write(current_iteration);
//
//         current_iteration++;
//     }
//     if (current_iteration == max_iterations)
//     {
//         throw std::runtime_error("Newton-Raphson iterations failed to converged.\n");
//     }
// }
}
