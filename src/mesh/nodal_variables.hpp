
#pragma once

#include "numeric/dense_matrix.hpp"

namespace neon
{
template <int variables>
class nodal_variables
{
public:
    static_assert(variables > 0, "Number of variables for the node must be greater than zero.");

    /// Storage type of the nodal variables (variables * number_of_nodes)
    using storage_type = matrixdx<variables>;

    static auto constexpr variables_per_node = variables;

public:
    /// Allocate storage for \p number_of_nodes
    explicit nodal_variables(std::int64_t const number_of_nodes)
        : u{storage_type::Zero(variables, number_of_nodes)}
    {
    }

    /// \return element variables based on \p local_indices
    template <typename indices_type>
    auto const view(indices_type const local_indices) const
    {
        return u(Eigen::placeholders::all, local_indices);
    }

    /// \return The total number of variables
    auto size() const noexcept { return u.size(); }

    /// Overwrite the nodal variables from a linear array
    void update(vector const& u_new) { vector::Map(u.data(), u.size()) = u_new; }

protected:
    /// Nodal variable storage
    storage_type u;
};
}
