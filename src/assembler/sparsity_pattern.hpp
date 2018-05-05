
#pragma once

#include "numeric/doublet.hpp"

#include <cstdint>
#include <vector>

namespace neon::fem
{
/// Compute the sparsity (nonzero) pattern of a sparse matrix and compress
/// the resulting data structure.  This function requires that the mesh_type
/// provide a \p local_dof_view method for each of the submeshes associated
/// with the \p mesh
/// This results in the non-zero entries in A set to zero.
template <typename sparse_matrix_type, typename mesh_type>
void compute_sparsity_pattern(sparse_matrix_type& A, mesh_type const& mesh)
{
    using integer_type = typename sparse_matrix_type::Index;

    static_assert(std::is_integral<integer_type>::value, "Index type must be an integer");

    std::vector<doublet<integer_type>> ij;
    ij.reserve(mesh.active_dofs());

    A.resize(mesh.active_dofs(), mesh.active_dofs());

    for (auto const& submesh : mesh.meshes())
    {
        // Loop over the elements and add in the non-zero components
        for (std::int64_t element{0}; element < submesh.elements(); element++)
        {
            auto const local_dof_view = submesh.local_dof_view(element);

            for (std::int64_t p{0}; p < local_dof_view.size(); p++)
            {
                for (std::int64_t q{0}; q < local_dof_view.size(); q++)
                {
                    ij.emplace_back(local_dof_view(p), local_dof_view(q));
                }
            }
        }
    }
    A.setFromTriplets(begin(ij), end(ij));
    A.finalize();
}
}
