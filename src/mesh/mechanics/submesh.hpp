
#pragma once

#include "mesh/basic_submesh.hpp"
#include "mesh/material_coordinates.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"

#include <memory>
#include <utility>
#include <tuple>

namespace neon::mechanics::detail
{
/// submesh is a Curiously Recurring Template Pattern class for enforcing
/// a minimum interface for all mechanics meshes at compile time.
template <class MeshType, class InternalVariableType>
class submesh : public basic_submesh
{
public:
    /// Alias to the fem mesh type
    using mesh_type = MeshType;

    /// Alias to the internal variable type
    using internal_variable_type = InternalVariableType;

public:
    explicit submesh(basic_submesh const& submesh) : basic_submesh(submesh) {}

    /// \return mapping of the element degrees of freedom to the process matrix
    auto local_dof_view(std::int64_t const element) const -> index_view
    {
        return static_cast<mesh_type*>(this)->local_dof_view(element);
    }

    /// Assemble the element stiffness matrix for a given \p element
    /// \param element The element number to assemble.
    /// \return the tangent consistent stiffness matrix
    auto tangent_stiffness(std::int64_t const element) const -> std::pair<index_view, matrix>
    {
        return static_cast<mesh_type*>(this)->tangent_stiffness(element);
    }

    /// \return the internal element force
    auto internal_force(std::int64_t const element) const -> std::pair<index_view, vector>
    {
        return static_cast<mesh_type*>(this)->internal_force(element);
    }

    /// \return the consistent mass matrix \sa diagonal_mass
    auto consistent_mass(std::int64_t const element) const -> std::pair<index_view, matrix>
    {
        return static_cast<mesh_type*>(this)->consistent_mass(element);
    }

    /// \return the diagonal mass matrix as a vector \sa consistent_mass
    auto diagonal_mass(std::int64_t const element) const -> std::pair<index_view, vector>
    {
        return static_cast<mesh_type*>(this)->diagonal_mass(element);
    }

    /// \return the number of degrees of freedom per node
    auto dofs_per_node() const noexcept -> std::int64_t
    {
        return static_cast<mesh_type*>(this)->dofs_per_node();
    }
};
}
