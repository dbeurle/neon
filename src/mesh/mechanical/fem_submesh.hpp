
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
/// fem_submesh is a Curiously Recurring Template Pattern class for enforcing
/// a minumum interface for all solid mechanics meshes at compile time.
template <class fem_mesh_type, class InternalVariableType>
class fem_submesh : public basic_submesh
{
public:
    using mesh_type = fem_mesh_type;
    using internal_variable_type = InternalVariableType;

public:
    explicit fem_submesh(basic_submesh const& submesh) : basic_submesh(submesh) {}

    /// \return mapping of the element degrees of freedom to the process matrix
    index_view local_dof_view(std::int64_t const element) const
    {
        return static_cast<mesh_type*>(this)->local_dof_view(element);
    }

    /// Assemble the element stiffness matrix for a given \p element
    /// \param element The element number to assemble.
    /// \return the tangent consistent stiffness matrix
    std::pair<index_view, matrix> tangent_stiffness(std::int64_t const element) const
    {
        return static_cast<mesh_type*>(this)->tangent_stiffness(element);
    }

    /// \return the internal element force
    std::pair<index_view, vector> internal_force(std::int64_t const element) const
    {
        return static_cast<mesh_type*>(this)->internal_force(element);
    }

    /// \return the consistent mass matrix \sa diagonal_mass
    std::pair<index_view, matrix> consistent_mass(std::int64_t const element) const
    {
        return static_cast<mesh_type*>(this)->consistent_mass(element);
    }

    /// \return the diagonal mass matrix as a vector \sa consistent_mass
    std::pair<index_view, vector> diagonal_mass(std::int64_t const element) const
    {
        return static_cast<mesh_type*>(this)->diagonal_mass(element);
    }

    /// \return the number of degrees of freedom per node
    auto dofs_per_node() const noexcept { return static_cast<mesh_type*>(this)->dofs_per_node(); }
};
}
