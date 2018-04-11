
#pragma once

#include <cstdint>
#include <vector>

namespace neon::fem
{
/// Apply homogeneous Dirichlet conditions to the matrix by zeroing the
/// row and column corresponding to the Dirichlet constrained degree of freedom
/// \param A system matrix
/// \param mesh The mesh containing boundary condition
template <typename sparse_matrix_type, typename mesh_type>
void apply_dirichlet_conditions(sparse_matrix_type& A, mesh_type const& mesh)
{
    // Keep track of entries we have visited as we walk to the
    // sparse matrix A
    std::vector<std::int32_t> non_zero_visitor;
    non_zero_visitor.reserve(A.nonZeros() / A.rows());

    for (auto const& [name, dirichlet_boundaries] : mesh.dirichlet_boundaries())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const fixed_dof : dirichlet_boundary.dof_view())
            {
                auto const diagonal_entry = A.coeff(fixed_dof, fixed_dof);

                non_zero_visitor.clear();

                for (typename sparse_matrix_type::InnerIterator it(A, fixed_dof); it; ++it)
                {
                    it.valueRef() = 0.0;

                    non_zero_visitor.push_back(A.IsRowMajor ? it.col() : it.row());
                }

                for (auto const& non_zero : non_zero_visitor)
                {
                    auto const row = A.IsRowMajor ? non_zero : fixed_dof;
                    auto const col = A.IsRowMajor ? fixed_dof : non_zero;

                    A.coeffRef(row, col) = 0.0;
                }
                // Reset the diagonal to the same value in an attempt to
                // preserve condition number
                A.coeffRef(fixed_dof, fixed_dof) = diagonal_entry;
            }
        }
    }
}

/// Apply homogeneous Dirichlet conditions to the matrix by zeroing the
/// row and column corresponding to the Dirichlet constrained degree of freedom
/// and correcting the right hand side vector for the contribution
/// \param A system matrix
/// \param x unknown vector
/// \param b right hand side
/// \param mesh The mesh containing boundary condition
template <typename sparse_matrix_type, typename vector_type, typename mesh_type>
void apply_dirichlet_conditions(sparse_matrix_type& A,
                                vector_type& x,
                                vector_type& b,
                                mesh_type const& mesh)
{
    // Keep track of entries we have visited as we walk to the
    // sparse matrix A
    std::vector<std::int32_t> non_zero_visitor;
    non_zero_visitor.reserve(A.nonZeros() / A.rows());

    for (auto const& [name, dirichlet_boundaries] : mesh.dirichlet_boundaries())
    {
        for (auto const& dirichlet_boundary : dirichlet_boundaries)
        {
            for (auto const fixed_dof : dirichlet_boundary.dof_view())
            {
                x(fixed_dof) = dirichlet_boundary.value_view();

                auto const diagonal_entry = A.coeff(fixed_dof, fixed_dof);

                non_zero_visitor.clear();

                for (typename sparse_matrix_type::InnerIterator it(A, fixed_dof); it; ++it)
                {
                    if (!A.IsRowMajor) b(it.row()) -= it.value() * x(fixed_dof);

                    it.valueRef() = 0.0;

                    non_zero_visitor.push_back(A.IsRowMajor ? it.col() : it.row());
                }

                for (auto const& non_zero : non_zero_visitor)
                {
                    auto const row = A.IsRowMajor ? non_zero : fixed_dof;
                    auto const col = A.IsRowMajor ? fixed_dof : non_zero;

                    if (A.IsRowMajor) b(row) -= A.coeff(row, col) * x(fixed_dof);

                    A.coeffRef(row, col) = 0.0;
                }
                // Reset the diagonal to the same value in an attempt to
                // preserve condition number
                A.coeffRef(fixed_dof, fixed_dof) = diagonal_entry;
                b(fixed_dof) = diagonal_entry * x(fixed_dof);
            }
        }
    }
}
}
