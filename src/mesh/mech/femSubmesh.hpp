
#pragma once

#include "mesh/Submesh.hpp"

#include "numeric/DenseMatrix.hpp"
#include "numeric/IndexTypes.hpp"

#include <tuple>

namespace neon::mech::detail
{
/**
 * femSubmesh is a Curiously Recurring Template Pattern class for enforcing
 * a minumum interface for all solid mechanics meshes at compile time.
 */
template <typename MeshType>
class femSubmesh : public Submesh
{
public:
    using mesh_type = MeshType;

public:
    femSubmesh<MeshType>(Submesh const& submesh) : Submesh(submesh) {}

    /** @return list of global degrees of freedom for an element */
    List const& local_dof_list(int64 const element) const
    {
        return static_cast<MeshType*>(this)->local_dof_list(element);
    }

    /** @return the tangent consistent stiffness matrix */
    [[nodiscard]] std::tuple<List const&, Matrix> tangent_stiffness(int64 const element) const
    {
        return static_cast<MeshType*>(this)->tangent_stiffness(element);
    }

    /** @return the internal element force */
    [[nodiscard]] std::tuple<List const&, Vector> internal_force(int64 const element) const
    {
        return static_cast<MeshType*>(this)->internal_force(element);
    }

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::tuple<List const&, Matrix> consistent_mass(int64 const element) const
    {
        return static_cast<MeshType*>(this)->consistent_mass(element);
    }

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::tuple<List const&, Vector> diagonal_mass(int64 const element) const
    {
        return static_cast<MeshType*>(this)->diagonal_mass(element);
    }

    [[nodiscard]] int64 dofs_per_node() const
    {
        return static_cast<MeshType*>(this)->dofs_per_node();
    }
};
}
