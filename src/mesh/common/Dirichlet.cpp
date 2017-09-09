
#include "Dirichlet.hpp"

namespace neon
{
Dirichlet::Dirichlet(List dofs, Json::Value const& times, Json::Value const& loads)
    : Boundary(times, loads), dofs(dofs)
{
}
}
