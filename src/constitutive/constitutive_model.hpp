
#pragma once

#include "internal_variables_forward.hpp"

#include <memory>
#include <set>

namespace neon
{
class material_property;

/// constitutive_model is the templated base class for all constitutive models.
/// The derived classes define their own internal variables, implement an update
/// internal variables routine and update a constitutive model for use in the global
/// assembly routine
template <int rank2_dimension, int rank4_dimension>
class constitutive_model
{
public:
    using internal_variable_t = internal_variables<rank2_dimension, rank4_dimension>;

public:
    /// Provide an internal variable class to be populated by the constitutive model
    explicit constitutive_model(std::shared_ptr<internal_variable_t>& variables)
        : variables(variables)
    {
    }

    virtual ~constitutive_model() = default;

    /// Update required internal variables and tangent matrix at quadrature points
    /// \param time_step_size Time step size (or load increment if quasi-static)
    virtual void update_internal_variables(double const time_step_size) = 0;

    /// \return A base class reference to the common material properties
    [[nodiscard]] virtual material_property const& intrinsic_material() const = 0;

    [[nodiscard]] virtual bool is_finite_deformation() const = 0;

    [[nodiscard]] virtual bool is_symmetric() const { return true; };

    [[nodiscard]] std::set<std::string> const& variable_names() const noexcept { return names; }

protected:
    /// Internal variables
    std::shared_ptr<internal_variable_t> variables;

    /// Internal variable names allocated by the constitutive model
    std::set<std::string> names;
};
}

#include "constitutive/constitutive_model_alias.hpp"
