
#include "dirichlet.hpp"

namespace neon::boundary
{
dirichlet::dirichlet(std::vector<int64> dofs, Json::Value const& times, Json::Value const& loads)
    : interpolator(times, loads), dofs(dofs)
{
}
}
