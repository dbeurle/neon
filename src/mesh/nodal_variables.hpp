
#pragma once

#include "numeric/dense_matrix.hpp"

namespace neon
{
template <int component_count>
class nodal_variables
{
public:
    static_assert(component_count > 0, "Number of variables for the node must be greater than zero.");

    /// Storage type of the nodal variables (component_count * number_of_nodes)
    using storage_type = matrixdx<component_count>;

    static auto constexpr components = component_count;

public:
    /// Allocate zeroed storage for \p number_of_nodes
    explicit nodal_variables(std::int64_t const number_of_nodes)
        : u{storage_type::Zero(components, number_of_nodes)}
    {
    }

    /// \return element variables based on \p local_indices
    template <typename indices_type>
    auto const view(indices_type const local_indices) const
    {
        return u(Eigen::placeholders::all, local_indices);
    }

    /// Obtain the underlying storage
    auto const& data() const noexcept { return u; }

    /// \return The total number of variables
    auto size() const noexcept { return u.cols(); }

    /// Overwrite the nodal variables from a linear array
    void update(vector const& u_new) { vector::Map(u.data(), u.size()) = u_new; }

    /// Overwrite the nodal variables from a two-dimensional array
    void update(storage_type const& u_new) { u = u_new; }

protected:
    /// Nodal variable storage
    storage_type u;
};
}
