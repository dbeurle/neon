
#include "nodal_value.hpp"

#include <utility>

namespace neon
{
nodal_value::nodal_value(std::vector<std::int32_t> dof_indices, json const& times, json const& loads)
    : boundary_condition{times, loads}, dof_indices{std::move(dof_indices)}
{
}

nodal_value::nodal_value(std::vector<std::int32_t> dof_indices,
                         json const& boundary_data,
                         std::string const& name,
                         double const generate_time_step)
    : boundary_condition{boundary_data, name, generate_time_step}, dof_indices{std::move(dof_indices)}
{
}

std::vector<std::int32_t> const& nodal_value::dof_view() const noexcept { return dof_indices; }

double nodal_value::value_view(double const load_factor) const
{
    return interpolate_prescribed_load(load_factor);
}
}
