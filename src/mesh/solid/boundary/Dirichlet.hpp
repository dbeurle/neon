
#pragma once

#include "numeric/DenseTypes.hpp"

#include <tuple>

namespace neon::solid
{
class Dirichlet
{
public:
    Dirichlet(List dofs, double const value, bool const is_ramped = true);

    auto const& dof_view() const { return dofs; }

    /** Get the value depending on the loading factor */
    auto const value_view(double const load_factor) const
    {
        return is_ramped ? (value_new - value_old) * load_factor + value_old : value_new;
    }

    /** Set the old value to the current value and update the current value */
    void update_value(double const value);

    /**
     * The boundary value is the same as the last step and are instantaneously applied
     */
    void inherit_from_last();

protected:
    List dofs;

    bool is_ramped;

    double value_old = 0.0;
    double value_new = 0.0;
};
}
