
#pragma once

#include "traits/mechanical.hpp"
#include "geometry/moment_inertia.hpp"

namespace neon::mechanical::beam
{
/// isotropic_linear is responsible for the calculation of the internal
/// variables for the linear beam theory.  This includes geometry of the beam
/// including the moment of inertia and the cross-sectional areas.
class isotropic_linear
{
public:
    using trait_type = mechanical::traits<theory::beam, discretisation::linear, true>;

public:
    isotropic_linear(json const& material_data, json const& section_data) {}

    void update_internal_variables() {}

    void update_stress() {}

protected:
    std::unique_ptr<geometry::section> section;

    isotropic_linear material;
};
}
