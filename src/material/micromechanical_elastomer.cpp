
#include "micromechanical_elastomer.hpp"

#include "numeric/float_compare.hpp"
#include "io/json.hpp"
#include "solver/time/euler_integration.hpp"
#include "solver/time/runge_kutta_integration.hpp"

#include <cmath>

namespace neon
{
micromechanical_elastomer::micromechanical_elastomer(json const& material_data)
    : isotropic_elastic_property(material_data)
{
    if (material_data.find("segments_per_chain") == end(material_data))
    {
        throw std::domain_error("segments_per_chain not specified in material data\n");
    }
    m_segments_per_chain = material_data["segments_per_chain"].get<double>();
}

ageing_micromechanical_elastomer::ageing_micromechanical_elastomer(json const& material_data)
    : micromechanical_elastomer(material_data)
{
    auto exception_string = [](std::string&& field) {
        return "\"" + field + "\" is not specified in \"material\" data";
    };

    if (material_data.find("scission_probability") == end(material_data))
    {
        throw std::domain_error(exception_string("scission_probability"));
    }
    if (material_data.find("recombination_probability") == end(material_data))
    {
        throw std::domain_error(exception_string("recombination_probability"));
    }
    if (material_data.find("cure_time") == end(material_data))
    {
        throw std::domain_error(exception_string("cure_time"));
    }
    if (material_data["cure_time"].get<double>() <= 10)
    {
        throw std::domain_error("The \"cure_time\" is less than 10.  This is likely a mistake.  "
                                "Please provide a value between 10 and 100");
    }

    scission = material_data["scission_probability"].get<double>();
    recombination = material_data["recombination_probability"].get<double>();
    m_shear_modulus_ia = (1.0 - material_data["cure_time"].get<double>() / 100.0) * shear_modulus();

    if (scission < 0.0 || recombination < 0.0)
    {
        throw std::domain_error("Material properties (probabilities) must be positive");
    }
}

double ageing_micromechanical_elastomer::creation_rate(double const active_shear_modulus,
                                                       double const inactive_shear_modulus,
                                                       double const active_segments,
                                                       double const inactive_segments) const
{
    // Inactive set recombination
    auto const alpha = 1.0 - std::pow(1.0 - recombination, inactive_segments + 1.0);
    // Active set combination
    auto const eta = 1.0 - std::pow(1.0 - recombination, active_segments + 1.0);

    return alpha * inactive_shear_modulus + 4 * eta * active_shear_modulus;
}

vector5 ageing_micromechanical_elastomer::integrate(vector5 const& z, double const time_step_size) const
{
    // Define the right hand side of the ageing evolution equations
    return runge_kutta_fourth_order(0.0, time_step_size, z, [=](double, vector5 y) -> vector5 {
        // Unpack the vector
        double const active_shear_modulus = y(0);
        double const inactive_shear_modulus = y(1);
        double const reduction = y(2);
        double const active_segments = y(3);
        double const inactive_segments = y(4);

        // Inactive set recombination
        auto const alpha = 1.0 - std::pow(1.0 - recombination, inactive_segments + 1.0);
        // Active set scission
        auto const beta = 1.0 - std::pow(1.0 - scission, active_segments);
        // Active set generation
        auto const eta = 1.0 - std::pow(1.0 - recombination, active_segments + 1.0);
        // Inactive set generation
        auto const nu = 1.0 - std::pow(1.0 - scission, inactive_segments);

        // Active set rate of change
        auto const y0 = alpha * inactive_shear_modulus + active_shear_modulus * (2.0 * eta - beta);
        // Inactive set rate of change
        auto const y1 = 2.0 * beta * active_shear_modulus
                        + inactive_shear_modulus * (nu - 2.0 * alpha);
        // Reduction factor rate of change
        auto const y2 = -reduction * (beta + 2.0 * eta);
        // Active set average segments per chain rate
        auto const y3 = (2.0 * alpha * inactive_shear_modulus * inactive_segments
                         - beta * active_segments * active_shear_modulus - active_segments * y0)
                        / active_shear_modulus;
        // Inactive set average segments per chain rate
        auto const y4 = inactive_shear_modulus > 1.0e-8
                            ? (beta * active_segments * active_shear_modulus
                               - 2.0 * alpha * inactive_segments * inactive_shear_modulus
                               - inactive_segments * y1)
                                  / inactive_shear_modulus
                            : 0.0;

        // Repack the data into y and return this as the update
        y(0) = y0;
        y(1) = y1;
        y(2) = y2;
        y(3) = y3;
        y(4) = y4;

        return y;
    });
}
}
