
#include "femDynamicMatrix.hpp"

#include "solver/linear/LinearSolver.hpp"

#include <chrono>
#include <json/value.h>
#include <termcolor/termcolor.hpp>

namespace neon::diffusion
{
femDynamicMatrix::femDynamicMatrix(femMesh& fem_mesh, Json::Value const& simulation_data)
    : femStaticMatrix(fem_mesh, simulation_data), time_solver(simulation_data["Time"])
{
    d = 250.0 * Vector::Ones(fem_mesh.active_dofs());
}

femDynamicMatrix::~femDynamicMatrix() { std::cout << "Dtor femdynamicmatrix" << std::endl; }

void femDynamicMatrix::solve()
{
    // Perform time dependent solution
    std::cout << "Solving " << fem_mesh.active_dofs() << " degrees of freedom\n";

    assemble_mass();

    assemble_stiffness();

    compute_external_force();

    while (time_solver.loop())
    {
        auto const start = std::chrono::high_resolution_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold
                  << "Time step: " << time_solver.iteration()
                  << ", time: " << time_solver.current_time() << termcolor::reset << std::endl;

        SparseMatrix A = M + time_solver.current_time_step_size() * K;

        Vector b = M * d + time_solver.current_time_step_size() * f;

        apply_dirichlet_conditions(A, d, b);

        linear_solver->solve(A, d, b);

        auto const end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Time step required " << elapsed_seconds.count() << "s\n";

        file_io.write(time_solver.iteration(), time_solver.current_time(), d);
    }
    std::cout << "Solver routine completed\n";
}

void femDynamicMatrix::assemble_mass()
{
    auto const start = std::chrono::high_resolution_clock::now();

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
    M.coeffs() = 0.0;

    doublets.clear();

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
}
