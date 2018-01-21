
#include "femDynamicMatrix.hpp"

#include "solver/linear/LinearSolver.hpp"

#include <chrono>
#include "io/json.hpp"
#include <termcolor/termcolor.hpp>

namespace neon::diffusion
{
femDynamicMatrix::femDynamicMatrix(femMesh& fem_mesh, json const& simulation_data)
    : femStaticMatrix(fem_mesh, simulation_data), time_solver(simulation_data["Time"])
{
    if (simulation_data.isMember("InitialConditions")
        && simulation_data["InitialConditions"].isMember("Uniform"))
    {
        d = simulation_data["InitialConditions"]["Uniform"].asDouble()
            * vector::Ones(fem_mesh.active_dofs());
    }
}

void femDynamicMatrix::solve()
{
    // Perform time dependent solution
    file_io.write(0, 0.0, d);

    std::cout << "\n"
              << std::string(4, ' ') << "Solving " << fem_mesh.active_dofs()
              << " degrees of freedom\n\n";

    assemble_mass();

    assemble_stiffness();

    compute_external_force();

    while (time_solver.loop())
    {
        auto const start = std::chrono::high_resolution_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold << "Time step "
                  << time_solver.iteration() << ", simulation time: " << time_solver.current_time()
                  << termcolor::reset << std::endl;

        SparseMatrix A = M + time_solver.current_time_step_size() * K;

        vector b = M * d + time_solver.current_time_step_size() * f;

        apply_dirichlet_conditions(A, d, b);

        linear_solver->solve(A, d, b);

        auto const end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Time step took " << elapsed_seconds.count() << "s\n";

        file_io.write(time_solver.iteration(), time_solver.current_time(), d);
    }
    std::cout << "Solver routine completed\n";
}

void femDynamicMatrix::assemble_mass()
{
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

    doublets.clear();

    auto const start = std::chrono::high_resolution_clock::now();

    for (auto const& submesh : fem_mesh.meshes())
    {
#pragma omp parallel for
        for (auto element = 0; element < submesh.elements(); ++element)
        {
            // auto const[dofs, m] = submesh.consistent_mass(element);
            auto const& tpl = submesh.consistent_mass(element);
            auto const& dofs = std::get<0>(tpl);
            auto const& m = std::get<1>(tpl);

            for (auto b = 0; b < dofs.size(); b++)
            {
                for (auto a = 0; a < dofs.size(); a++)
                {
#pragma omp atomic
                    M.coeffRef(dofs[a], dofs[b]) += m(a, b);
                }
            }
        }
    }

    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Mass assembly took " << elapsed_seconds.count() << "s\n";
}
}
