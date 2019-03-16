
#pragma once

/// @file

#include "constitutive/internal_variables.hpp"
#include "geometry/profile.hpp"
#include "io/json_forward.hpp"
#include "material/isotropic_elastic_property.hpp"
#include "traits/mechanics.hpp"

namespace neon::mechanics::beam
{
/// isotropic_linear is responsible for the calculation of the internal
/// variables for the linear beam theory.  This includes geometry of the beam
/// including the moment of inertia and the cross-sectional areas.
class isotropic_linear
{
public:
    /// Trait for beam
    using traits = mechanics::traits<theory::beam, discretisation::linear, true>;

    /// Type alias for internal variables
    using variable_type = internal_variables<traits::rank_two_tensor, traits::rank_four_tensor>;

public:
    isotropic_linear(std::shared_ptr<variable_type>& variables, json const& material_data);

    void update_internal_variables();

    void update_stress();

protected:
    std::shared_ptr<variable_type> variables;

    isotropic_elastic_property material;
};
}
