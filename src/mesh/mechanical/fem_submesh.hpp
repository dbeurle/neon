
#pragma once

#include "mesh/basic_submesh.hpp"
#include "mesh/material_coordinates.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"

#include <memory>
#include <utility>
#include <tuple>

namespace neon::mechanical::detail
{
/**
 * fem_submesh is a Curiously Recurring Template Pattern class for enforcing
 * a minumum interface for all solid mechanics meshes at compile time.
 */
template <class MeshType, class InternalVariableType>
class fem_submesh : public basic_submesh
{
public:
    using mesh_type = MeshType;
    using internal_variable_type = InternalVariableType;

public:
    explicit fem_submesh(basic_submesh const& submesh) : basic_submesh(submesh) {}

    /** @return mapping of the element degrees of freedom to the process matrix */
    local_indices const& local_dof_view(std::int32_t const element) const
    {
        return static_cast<mesh_type*>(this)->local_dof_view(element);
    }

    /** @return the tangent consistent stiffness matrix */
    [[nodiscard]] std::pair<local_indices const&, matrix> tangent_stiffness(
        std::int32_t const element) const {
        return static_cast<mesh_type*>(this)->tangent_stiffness(element);
    }

    /** @return the internal element force */
    std::pair<local_indices const&, vector> internal_force(std::int32_t const element) const
    {
        return static_cast<mesh_type*>(this)->internal_force(element);
    }

    /** @return the consistent mass matrix \sa diagonal_mass */
    std::pair<local_indices const&, matrix> consistent_mass(std::int32_t const element) const
    {
        return static_cast<mesh_type*>(this)->consistent_mass(element);
    }

    /** @return the diagonal mass matrix as a vector \sa consistent_mass */
    std::pair<local_indices const&, vector> diagonal_mass(std::int32_t const element) const
    {
        return static_cast<mesh_type*>(this)->diagonal_mass(element);
    }

    /** @return the number of degrees of freedom per node */
    std::int32_t dofs_per_node() const { return static_cast<MeshType*>(this)->dofs_per_node(); }

protected:
};
}
