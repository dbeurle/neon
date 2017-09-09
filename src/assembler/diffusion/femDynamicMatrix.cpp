//
// #include "femDynamicMatrix.hpp"
//
// #include "solver/linear/LinearSolver.hpp"
//
// #include <chrono>
// #include <termcolor/termcolor.hpp>
//
// namespace neon::diffusion
// {
// femDynamicMatrix::femDynamicMatrix(femMesh& fem_mesh,
//                                    Visualisation&& visualisation,
//                                    Json::Value const& solver_data,
//                                    Json::Value const& time_data)
//     : femStaticMatrix(fem_mesh, std::forward<Visualisation>(visualisation), solver_data,
//     time_data),
//       a(Vector::Zero(fem_mesh.active_dofs())),
//       v(Vector::Zero(fem_mesh.active_dofs())),
//       newmark(time_data)
// {
// }
//
// void femDynamicMatrix::solve()
// {
//     // Perform Newton-Raphson iterations
//     std::cout << "Solving " << fem_mesh.active_dofs() << " degrees of freedom\n";
//
//     assemble_mass();
//
//     apply_displacement_boundaries();
//
//     // TODO Put the internal variable update here
//
//     std::cout << "Initial displacements\n" << d << std::endl;
//
//     compute_internal_force();
//     std::cout << "Initial internal force vector\n" << fint << std::endl;
//
//     auto f = -fint;
//
//     // Compute initial accelerations (a = inv(M) * f) where f = fext - fint
//     a = M.cwiseInverse().cwiseProduct(f);
//
//     std::cout << "Initial trial accelerations\n" << a << std::endl;
//
//     // Estimate solution
//     // d = newmark.approximate_displacements(a, v, d);
//
//     std::cout << "Newmark trial displacements\n" << d << std::endl;
//
//     while (newmark.time_loop())
//     {
//         std::cout << "Performing time step\n";
//
//         perform_equilibrium_iterations();
//
//         // fem_mesh.write();
//     }
//     std::cout << "Solver routine finished\n";
// }
//
// void femDynamicMatrix::assemble_mass()
// {
//     auto start = std::chrono::high_resolution_clock::now();
//
//     M = Vector::Zero(fem_mesh.active_dofs());
//
//     for (auto const& submesh : fem_mesh.meshes())
//     {
//         for (auto element = 0; element < submesh.elements(); ++element)
//         {
//             auto const[dofs, m] = submesh.diagonal_mass(element);
//
//             for (auto a = 0; a < dofs.size(); a++)
//             {
//                 M(dofs[a]) += m(a);
//             }
//         }
//     }
//
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//
//     std::cout << "  Assembly of mass matrix took " << elapsed_seconds.count() << "s\n";
// }
//
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
// }
