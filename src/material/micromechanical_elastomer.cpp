
#include "micromechanical_elastomer.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/float_compare.hpp"

#include "io/json.hpp"

#include <cmath>

using namespace neon;

micromechanical_elastomer::micromechanical_elastomer(json const& material_data)
    : isotropic_elastic_property(material_data)
{
    if (!material_data.count("SegmentsPerChain"))
    {
        throw std::domain_error("SegmentsPerChain not specified in material data\n");
    }
}

ageing_micromechanical_elastomer::ageing_micromechanical_elastomer(json const& material_data)
    : micromechanical_elastomer(material_data)
{
    auto exception_string = [](std::string&& field) {
        return "\"" + field + "\" is not specified in \"Material\" data";
    };

    if (!material_data.count("SegmentDecayRate"))
    {
        throw std::domain_error(exception_string("SegmentDecayRate"));
    }
    if (!material_data.count("ScissionProbability"))
    {
        throw std::domain_error(exception_string("ScissionProbability"));
    }
    if (!material_data.count("CrosslinkGrowthRate"))
    {
        throw std::domain_error(exception_string("CrosslinkGrowthRate"));
    }

    segment_decay_rate = material_data["SegmentDecayRate"];
    scission_probability = material_data["ScissionProbability"];
    crosslink_growth_rate = material_data["CrosslinkGrowthRate"];

    if (scission_probability < 0.0 || crosslink_growth_rate < 0.0 || segment_decay_rate < 0.0)
    {
        throw std::domain_error("Material properties for the segments must be positive");
    }
}

double ageing_micromechanical_elastomer::compute_new_segment(double const current_segment,
                                                             double const time_step_size) const
{
    return current_segment * (1.0 - segment_decay_rate * time_step_size);
}

double ageing_micromechanical_elastomer::compute_new_shear_modulus(double const time_step_size) const
{
    return crosslink_growth_rate * time_step_size;
}

std::vector<double> ageing_micromechanical_elastomer::compute_shear_moduli(
    std::vector<double> current_shear_moduli,
    std::vector<double> const& current_segments,
    double const time_step_size) const
{
    std::transform(std::begin(current_shear_moduli),
                   std::end(current_shear_moduli),
                   std::begin(current_segments),
                   std::begin(current_shear_moduli),
                   [this, time_step_size](auto const G, auto const N) {
                       return G * (1.0 - std::pow(1.0 - scission_probability, N) * time_step_size);
                   });

    return current_shear_moduli;
}
