
#pragma once

/// @file

#include "mesh/boundary/boundary.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
/// nodal_value handles a prescribed nodal value for potentially multiple
/// degrees of freedom with the same value history prescribed.  This can be
/// used to represent nodal based boundary conditions
/// \sa dirichlet
class nodal_value : public boundary
{
public:
    explicit nodal_value(std::vector<std::int32_t> dof_indices, json const& times, json const& loads);

    explicit nodal_value(std::vector<std::int32_t> dof_indices,
                         json const& boundary_data,
                         std::string const& name,
                         double const generate_time_step);

    /// \return view of the unique dof indices for this boundary
    [[nodiscard]] std::vector<std::int32_t> const& dof_view() const noexcept;

    /// \return boundary value depending on the loading factor
    [[nodiscard]] double value_view(double const load_factor = 1.0) const;

protected:
    /// Unique dof indices to apply the value
    std::vector<std::int32_t> dof_indices;
};
}
