
#pragma once

#include <json/forwards.h>

#include "numeric/DenseTypes.hpp"
#include "numeric/SparseTypes.hpp"
#include "solver/TimeControl.hpp"

namespace neon
{
class NewmarkBeta
{
public:
    NewmarkBeta(Json::Value const& time_solver_data);

    double mass_scaling_factor() const;

    Vector approximate_displacements(Vector const& a, Vector const& v, Vector const& d) const;

    Vector approximate_velocities(Vector const& a, Vector const& v) const;

    Vector accelerations(Vector const& a, Vector const& v, Vector const& d) const;

    Vector velocities(Vector const& a, Vector const& v) const;

    /** Perform the time integration until returns false */
    bool time_loop();

protected:
    bool are_parameters_unstable() const;

protected:
    TimeControl time_control;
    double artifical_viscosity;
    double beta_parameter;
};

inline double NewmarkBeta::mass_scaling_factor() const
{
    auto const Δt = time_control.current_time_step_size();
    auto const β = beta_parameter;
    return 1.0 / (β * std::pow(Δt, 2));
}

inline Vector NewmarkBeta::approximate_displacements(Vector const& a,
                                                     Vector const& v,
                                                     Vector const& d) const
{
    auto const Δt = time_control.current_time_step_size();
    auto const β = beta_parameter;

    return d + Δt * v + std::pow(Δt, 2) / 2.0 * (1.0 - 2.0 * β) * a;
}

inline Vector NewmarkBeta::approximate_velocities(Vector const& a, Vector const& v) const
{
    auto const Δt = time_control.current_time_step_size();
    auto const γ = artifical_viscosity;

    return v + (1.0 - γ) * Δt * a;
}

inline Vector NewmarkBeta::accelerations(Vector const& a, Vector const& v, Vector const& d) const
{
    auto const Δt = time_control.current_time_step_size();
    auto const β = beta_parameter;

    return 1.0 / (β * std::pow(Δt, 2)) * (d - approximate_displacements(a, v, d));
}
inline Vector NewmarkBeta::velocities(Vector const& a, Vector const& v) const
{
    auto const Δt = time_control.current_time_step_size();
    auto const γ = artifical_viscosity;

    return approximate_velocities(a, v) + γ * Δt * a;
}
}
