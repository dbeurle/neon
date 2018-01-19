
#pragma once

#include "mesh/Submesh.hpp"

#include "mesh/MaterialCoordinates.hpp"

#include "numeric/DenseMatrix.hpp"
#include "numeric/IndexTypes.hpp"

#include <memory>
#include <tuple>

namespace neon::mechanical::detail
{
/**
 * femSubmesh is a Curiously Recurring Template Pattern class for enforcing
 * a minumum interface for all solid mechanics meshes at compile time.
 */
template <class MeshType, class InternalVariableType>
class femSubmesh : public Submesh
{
public:
    using mesh_type = MeshType;
    using internal_variable_type = InternalVariableType;

public:
    explicit femSubmesh(Submesh const& submesh) : Submesh(submesh) {}

    /** @return mapping of the element degrees of freedom to the process matrix */
    local_indices const& local_dof_view(int32 const element) const
    {
        return static_cast<mesh_type*>(this)->local_dof_view(element);
    }

    /** @return the tangent consistent stiffness matrix */
    [[nodiscard]] std::pair<local_indices const&, matrix> tangent_stiffness(int32 const element) const
    {
        return static_cast<mesh_type*>(this)->tangent_stiffness(element);
    }

    /** @return the internal element force */
    std::pair<local_indices const&, vector> internal_force(int32 const element) const
    {
        return static_cast<mesh_type*>(this)->internal_force(element);
    }

    /** @return the consistent mass matrix \sa diagonal_mass */
    std::pair<local_indices const&, matrix> consistent_mass(int32 const element) const
    {
        return static_cast<mesh_type*>(this)->consistent_mass(element);
    }

    /** @return the diagonal mass matrix as a vector \sa consistent_mass */
    std::pair<local_indices const&, vector> diagonal_mass(int32 const element) const
    {
        return static_cast<mesh_type*>(this)->diagonal_mass(element);
    }

    /** @return the number of degrees of freedom per node */
    int32 dofs_per_node() const { return static_cast<MeshType*>(this)->dofs_per_node(); }

protected:
};
}
