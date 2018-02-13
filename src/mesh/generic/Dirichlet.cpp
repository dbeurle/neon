
#include "Dirichlet.hpp"

namespace neon
{
Dirichlet::Dirichlet(local_indices dofs, json const& times, json const& loads)
    : Boundary(times, loads), dofs(dofs)
{
}

Dirichlet::Dirichlet(local_indices dofs,
                     json const& boundary,
                     std::string const& name,
                     double const generate_time_step)
    : Boundary(boundary, name, generate_time_step), dofs(dofs)
{
}
}
