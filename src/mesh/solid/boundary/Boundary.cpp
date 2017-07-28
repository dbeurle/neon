
#include "Boundary.hpp"

namespace neon
{
Boundary::Boundary(double const prescribed_value, bool const is_load_ramped)
    : is_load_ramped(is_load_ramped), value_new(prescribed_value)
{
}

void Boundary::internal_restart(double const prescribed_value_new, bool const is_load_ramped)
{
    value_old = value_new;
    value_new = prescribed_value_new;
    this->is_load_ramped = is_load_ramped;
}

void Boundary::internal_restart()
{
    value_old = value_new;
    is_load_ramped = false;
}
}
