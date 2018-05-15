
#include "nodal_value.hpp"

namespace neon
{
nodal_value::nodal_value(std::vector<std::int32_t> unique_dofs, json const& times, json const& loads)
    : boundary{times, loads}, unique_dofs{unique_dofs}
{
}

nodal_value::nodal_value(std::vector<std::int32_t> unique_dofs,
                         json const& boundary_data,
                         std::string const& name,
                         double const generate_time_step)
    : boundary{boundary_data, name, generate_time_step}, unique_dofs{unique_dofs}
{
}
}
