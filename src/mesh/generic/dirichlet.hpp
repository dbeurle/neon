
#pragma once

#include "mesh/generic/interpolator.hpp"

#include "numeric/IndexTypes.hpp"

namespace neon::boundary
{
class dirichlet : public interpolator
{
public:
    explicit dirichlet(std::vector<int64> dofs, Json::Value const& times, Json::Value const& loads);

    [[nodiscard]] auto const& dof_view() const { return dofs; }

    /** Get the value depending on the loading factor */
    [[nodiscard]] auto value_view(double const load_factor = 1.0) const
    {
        return interpolate_prescribed_load(load_factor);
    }

protected:
    std::vector<int64> dofs;
};
}
