
#pragma once

#include "math/view.hpp"
#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"

#include <cstdint>
#include <vector>
#include <utility>

namespace neon
{
/// recovery is the base class for all element recovery methods for C0 finite
/// elements.  This defines the interface to handle the extrapolation of the
/// derivative values (dual variables) to the nodes for visualisation.  Direct
/// visualisation of the integration point values is not possible with the current
/// VTK visualisation code and will require projection to the nodes until support
/// is added in VTK.
class recovery
{
public:
    using value_type = std::pair<vector, std::vector<std::uint16_t>>;

public:
    virtual ~recovery() = default;

    /// Perform the projection on the flattened array of variables
    /// \return a pair of the values and the insertion count
    /// \param variables A flat array of scalars
    /// \param stride The encoding along for each group of integration points
    /// \param node_count Total size of each node group
    virtual auto project(std::vector<double> const& variables,
                         stride_view<> const& stride,
                         indices const& node_indices,
                         std::int64_t const node_count) const noexcept -> value_type = 0;
};

/// Compute the extrapolation matrix to allow for quadrature valued variables
/// to be averaged to the nodal points without ill-effects when using a
/// least squares (for example with quadratric tetrahedron elements)
/// using the algorithm in \cite Durand2014
class local_extrapolation : public recovery
{
public:
    local_extrapolation(matrix const& shape_functions,
                        matrix const& local_nodal_coordinates,
                        matrix const& local_quadrature_coordinates);

    auto project(std::vector<double> const& variables,
                 stride_view<> const& stride,
                 indices const& node_indices,
                 std::int64_t const node_count) const noexcept -> value_type override final;

    auto extrapolation_matrix() const noexcept -> matrix const& { return m_extrapolation; }

private:
    void compute_extrapolation_matrix(matrix const& shape_functions,
                                      matrix const& local_nodal_coordinates,
                                      matrix const& local_quadrature_coordinates);

private:
    /// Quadrature point to nodal point mapping
    matrix m_extrapolation;
};

namespace internal
{
auto linear_basis(double const x, double const y) -> vector3 { return vector3(1.0, x, y); }

auto linear_basis(double const x, double const y, double const z) -> vector4
{
    return vector4(1.0, x, y, z);
}
}

/// superconvergent_recovery is the classical superconvergent patch recovery
/// method by Zienkiewicz and Zhu.  Originally designed for gradient recovery
/// (stress) values at quadrature points by fitting a characteristic polynomial
/// of a patch around an element
template <typename MeshType>
auto superconvergent_patch_recovery(std::int64_t const node_index, MeshType const& mesh) -> double
{
    // Obtain the highest basis function
    auto basis = internal::linear_basis();

    // for each node check if it exists in a submesh

    // Compute the coordinates at the quadrature point
    // coordinates * shape_functions

    return 0.0;
}
}
