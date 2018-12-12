
#pragma once

#include "numeric/sparse_matrix.hpp"

#include "assembler/sparsity_pattern.hpp"
#include "assembler/homogeneous_dirichlet.hpp"
#include "solver/eigen/arpack.hpp"

#include <chrono>
#include <iostream>

namespace neon::mechanics
{
/// linear_buckling_matrix assembles and solves the eigenvalue buckling problem
/// for linear constitutive models only.
template <typename MeshType>
class linear_buckling_matrix
{
public:
    using mesh_type = MeshType;

public:
    linear_buckling_matrix(mesh_type& mesh)
        : mesh(mesh), solver(5, eigen_solver::eigen_spectrum::lower)
    {
    }

    /// Compute the eigenvalue for the buckling load and the corresponding
    /// buckling mode.
    void solve();

protected:
    void assemble_stiffness();

protected:
    /// fem mesh
    mesh_type& mesh;
    /// Stiffness matrix
    sparse_matrix K;
    /// Eigenvalue solver
    arpack solver;
};

template <typename MeshType>
void linear_buckling_matrix<MeshType>::solve()
{
    // Solve the eigenvalue problem for the first eigenvalue only
    assemble_stiffness();

    // apply_dirichlet_conditions(K, mesh);

    // solver.solve(K);

    // file_io.write(0, 0.0, d);
}

template <typename MeshType>
void linear_buckling_matrix<MeshType>::assemble_stiffness()
{
    compute_sparsity_pattern(K, mesh);

    auto const start = std::chrono::steady_clock::now();

    K.coeffs() = 0.0;

    for (auto const& submesh : mesh.meshes())
    {
        for (std::int64_t element = 0; element < submesh.elements(); ++element)
        {
            auto const [dofs, ke] = submesh.tangent_stiffness(element);

            for (std::int64_t a{0}; a < dofs.size(); a++)
            {
                for (std::int64_t b{0}; b < dofs.size(); b++)
                {
                    K.coefficient_update(dofs(a), dofs(b), ke(a, b));
                }
            }
        }
    }

    auto const end = std::chrono::steady_clock::now();

    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Stiffness assembly took " << elapsed_seconds.count()
              << "s\n";
}
}
