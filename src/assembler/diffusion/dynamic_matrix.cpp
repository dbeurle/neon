
#include "dynamic_matrix.hpp"

#include "assembler/homogeneous_dirichlet.hpp"
#include "solver/linear/linear_solver.hpp"
#include "numeric/doublet.hpp"
#include "io/json.hpp"

#include <termcolor/termcolor.hpp>
#include <tbb/parallel_for.h>

#include <chrono>

namespace neon::diffusion
{
dynamic_matrix::dynamic_matrix(mesh_type& mesh, json const& simulation_data)
    : static_matrix(mesh, simulation_data), time_solver(simulation_data["time"])
{
    if (simulation_data.find("initial_conditions") != simulation_data.end()
        && simulation_data["initial_conditions"].find("uniform")
               != simulation_data["initial_conditions"].end())
    {
        d = simulation_data["initial_conditions"]["uniform"].get<double>()
            * vector::Ones(mesh.active_dofs());
    }
    else
    {
        // throw std::domain_error(
        //     "\"initial_conditions\" must be provided for a transient simulation.");
    }
}

void dynamic_matrix::solve()
{
    // Perform time dependent solution
    // file_io.write(0, 0.0, d);

    std::cout << "\n"
              << std::string(4, ' ') << "Solving " << mesh.active_dofs()
              << " degrees of freedom\n\n";

    assemble_mass();

    assemble_stiffness();

    compute_external_force();

    sparse_matrix A;
    vector b;

    while (time_solver.loop())
    {
        auto const start = std::chrono::steady_clock::now();

        std::cout << std::string(4, ' ') << termcolor::blue << termcolor::bold << "Time step "
                  << time_solver.iteration() << ", simulation time: " << time_solver.current_time()
                  << termcolor::reset << std::endl;

        A = M + time_solver.current_time_step_size() * K;

        b = M * d + time_solver.current_time_step_size() * f;

        apply_dirichlet_conditions(A, d, b, mesh);

        solver->solve(A, d, b);

        auto const end = std::chrono::steady_clock::now();
        std::chrono::duration<double> const elapsed_seconds = end - start;
        std::cout << std::string(6, ' ') << "Time step took " << elapsed_seconds.count() << "s\n";

        // mesh.write(time_solver.iteration(), time_solver.current_time());
    }
    std::cout << "Solver routine completed\n";
}

void dynamic_matrix::assemble_mass()
{
    M.resize(mesh.active_dofs(), mesh.active_dofs());

    std::vector<doublet<int>> doublets;
    doublets.reserve(mesh.active_dofs());

    for (auto const& submesh : mesh.meshes())
    {
        for (std::int64_t element = 0; element < submesh.elements(); element++)
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
    M.setFromTriplets(doublets.begin(), doublets.end());

    doublets.clear();

    auto const start = std::chrono::steady_clock::now();

    for (auto const& submesh : mesh.meshes())
    {
        for (std::int64_t element = 0; element < submesh.elements(); ++element)
        {
            auto const& [dofs, m] = submesh.consistent_mass(element);

            for (std::int64_t a{0}; a < dofs.size(); a++)
            {
                for (std::int64_t b{0}; b < dofs.size(); b++)
                {
                    M.add_to(dofs(a), dofs(b), m(a, b));
                }
            }
        }
    }

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Mass assembly took " << elapsed_seconds.count() << "s\n";
}
}
