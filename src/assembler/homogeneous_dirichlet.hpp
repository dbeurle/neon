
#pragma once

#include <cstdint>
#include <vector>

namespace neon::fem
{
/// Apply dirichlet conditions to the system defined by A, x, and b.
/// This method selects the row and column of the degree of freedom with an
/// imposed Dirichlet condition.
/// The algorithm takes the row and performs a zeroing operation on it.  To
/// satisfy the equation system, the constrained column is multiplied with
/// the Dirichlet value and subtracted from the matching DoF on the right
/// hand side vector.  The column is then zeroed and the diagonal DoF is
/// then corrected such that \f$ A_{dof} * x_{dof} = f_{dof} == A_{dof} * x_{dof} \f$ so the
/// equation system is satisfied.
/// For inner and outer vector reference see
/// https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
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
