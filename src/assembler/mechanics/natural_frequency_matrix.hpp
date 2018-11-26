
#pragma once

#include "numeric/sparse_matrix.hpp"

#include "assembler/sparsity_pattern.hpp"
#include "assembler/homogeneous_dirichlet.hpp"
#include "solver/eigen/eigenvalue_solver.hpp"

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
    natural_frequency_matrix(mesh_type&& mesh) : mesh(std::move(mesh)), solver(5) {}

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
    eigenvalue_solver solver;

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

    solver.solve(K, M);

    print_eigenvalue_table();

    mesh.write(solver.eigenvalues().cwiseSqrt(), solver.eigenvectors());
}

template <typename MeshType>
void natural_frequency_matrix<MeshType>::assemble_stiffness()
{
    fem::compute_sparsity_pattern(K, mesh);

    auto const start = std::chrono::steady_clock::now();

    K.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
        for (std::int64_t element = 0; element < submesh.elements(); ++element)
        {
            auto const& [dofs, local_tangent] = submesh.tangent_stiffness(element);

            for (std::int64_t a{0}; a < dofs.size(); a++)
            {
                for (std::int64_t b{0}; b < dofs.size(); b++)
                {
                    K.coefficient_update(dofs(a), dofs(b), local_tangent(a, b));
                }
            }
        }
    }

    auto const end = std::chrono::steady_clock::now();

    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Stiffness matrix assembly took " << elapsed_seconds.count()
              << "s\n";
}

template <typename MeshType>
void natural_frequency_matrix<MeshType>::assemble_mass()
{
    fem::compute_sparsity_pattern(M, mesh);

    auto const start = std::chrono::steady_clock::now();

    M.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
        for (std::int64_t element = 0; element < submesh.elements(); ++element)
        {
            auto const& [dofs, local_mass] = submesh.consistent_mass(element);

            for (std::int64_t a{0}; a < dofs.size(); a++)
            {
                for (std::int64_t b{0}; b < dofs.size(); b++)
                {
                    M.coefficient_update(dofs(a), dofs(b), local_mass(a, b));
                }
            }
        }
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

    for (std::int64_t index{0}; index < solver.eigenvalues().size(); ++index)
    {
        std::cout << std::string(indent, ' ') << std::setw(10) << index + 1
                  << std::setw(17)
                  // print radians / time
                  << std::sqrt(solver.eigenvalues()(index))
                  << std::setw(16)
                  // print cycles / time
                  << std::sqrt(solver.eigenvalues()(index)) / (2.0 * M_PI) << "\n";
    }
}
}
