
#pragma once

/// @file

#include "numeric/sparse_matrix.hpp"

#include "assembler/sparsity_pattern.hpp"
#include "assembler/homogeneous_dirichlet.hpp"
#include "solver/eigen/arpack.hpp"

#include <chrono>
#include <iostream>

namespace neon::mechanics
{
/// natural_frequency_matrix assembles and solves the eigenvalue buckling problem
/// for linear constitutive models only.
template <typename MeshType>
class natural_frequency_matrix
{
public:
    using mesh_type = MeshType;

public:
    natural_frequency_matrix(mesh_type&& mesh, std::unique_ptr<eigen_solver>&& solver)
        : mesh(std::move(mesh)), solver{std::move(solver)}
    {
    }

    /// Compute the eigenvalues corresponding to the resonant frequency of the
    /// structure.
    void solve();

protected:
    void assemble_stiffness();

    void assemble_mass();

protected:
    /// fem mesh
    mesh_type mesh;

    /// Stiffness matrix
    sparse_matrix K;
    /// Mass matrix
    sparse_matrix M;

    /// Eigenvalue solver
    std::unique_ptr<eigen_solver> solver;

private:
    void print_eigenvalue_table() const;
};

template <typename MeshType>
void natural_frequency_matrix<MeshType>::solve()
{
    assemble_stiffness();

    assemble_mass();

    apply_dirichlet_conditions(K, mesh);
    apply_dirichlet_conditions(M, mesh);

    solver->solve(K, M);

    print_eigenvalue_table();

    mesh.write(solver->eigenvalues().cwiseSqrt(), solver->eigenvectors());
}

template <typename MeshType>
void natural_frequency_matrix<MeshType>::assemble_stiffness()
{
    compute_sparsity_pattern(K, mesh);

    auto const start = std::chrono::steady_clock::now();

    K.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
        tbb::parallel_for(std::int64_t{0}, submesh.elements(), [&](auto const element) {
            auto const dof_indices = submesh.local_dof_view(element);
            auto const& local_tangent = submesh.tangent_stiffness(element);

            for (std::int64_t a{0}; a < dof_indices.size(); a++)
            {
                for (std::int64_t b{0}; b < dof_indices.size(); b++)
                {
                    K.add_to(dof_indices(a), dof_indices(b), local_tangent(a, b));
                }
            }
        });
    }

    auto const end = std::chrono::steady_clock::now();

    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Stiffness matrix assembly took " << elapsed_seconds.count()
              << "s\n";
}

template <typename MeshType>
void natural_frequency_matrix<MeshType>::assemble_mass()
{
    compute_sparsity_pattern(M, mesh);

    auto const start = std::chrono::steady_clock::now();

    M.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
        tbb::parallel_for(std::int64_t{0}, submesh.elements(), [&](auto const element) {
            auto const dof_indices = submesh.local_dof_view(element);
            auto const& local_mass = submesh.consistent_mass(element);

            for (std::int64_t a{0}; a < dof_indices.size(); a++)
            {
                for (std::int64_t b{0}; b < dof_indices.size(); b++)
                {
                    M.add_to(dof_indices(a), dof_indices(b), local_mass(a, b));
                }
            }
        });
    }

    auto const end = std::chrono::steady_clock::now();

    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Mass matrix assembly took " << elapsed_seconds.count()
              << "s\n";
}

template <typename MeshType>
void natural_frequency_matrix<MeshType>::print_eigenvalue_table() const
{
    constexpr auto indent = 6;
    std::cout << '\n'
              << std::string(indent, ' ') << "Eigenvalue | "
              << "radians / time | "
              << "cycles / time\n"
              << std::string(indent, ' ') << std::setfill('-') << std::setw(44) << '\n'
              << std::setfill(' ');

    for (std::int64_t index{0}; index < solver->eigenvalues().size(); ++index)
    {
        std::cout << std::string(indent, ' ') << std::setw(10) << index + 1
                  << std::setw(17)
                  // print radians / time
                  << std::sqrt(solver->eigenvalues()(index))
                  << std::setw(16)
                  // print cycles / time
                  << std::sqrt(solver->eigenvalues()(index)) / (2.0 * M_PI) << "\n";
    }
}
}
