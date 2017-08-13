
#include "MicromechanicalElastomer.hpp"

#include "numeric/DenseTypes.hpp"

#include <json/value.h>

#include <cmath>
#include <iostream>
#include <random>

#include <range/v3/action.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace neon
{
/**
 * Computes the probability mass function of the Binomial distribution
 * @param n - The number of trials
 * @param k - Number of successes
 * @param p - Probability of success
 * @return The probabilty Pr(k;n,p)
 */
double binomial_pmf(int const n, int const k, double const p = 0.5)
{
    return std::exp(std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1))
           * std::pow(p, k) * std::pow(1.0 - p, n - k);
}

MicromechanicalElastomer::MicromechanicalElastomer(Json::Value const& material_data)
    : LinearElastic(material_data)
{
    if (!material_data.isMember("Segments"))
    {
        throw std::runtime_error("Segments not specified in material data\n");
    }
    compute_chains_and_segments(material_data["Segments"]);
}

void MicromechanicalElastomer::compute_chains_and_segments(Json::Value const& segments_data)
{
    // Basic error checking
    if (!segments_data.isMember("Groups"))
    {
        throw std::runtime_error("Groups not specified in \"Segment\" data\n");
    }
    if (!segments_data.isMember("Average"))
    {
        throw std::runtime_error("Average not specified in \"Segment\" data\n");
    }
    if (!segments_data.isMember("StandardDeviation"))
    {
        throw std::runtime_error("StandardDeviation not specified in \"Segment\" data\n");
    }
    if (!segments_data.isMember("ScissionLikelihood"))
    {
        throw std::runtime_error("ScissionLikelihood not specified in \"Segment\" "
                                 "data\n");
    }

    number_of_groups = segments_data["Groups"].asInt();
    p_scission = segments_data["ScissionLikelihood"].asDouble();

    auto const N_avg = segments_data["Average"].asInt();
    auto const N_std = segments_data["StandardDeviation"].asInt();

    auto const n0 = LinearElastic::shear_modulus() / (boltzmann_constant * temperature);

    // Create a normal distribution generator and sample the same number
    // with the segments per chain as the mean value
    std::random_device random_dev;
    std::mt19937 generator(random_dev());

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> normal_dist(N_avg, N_std);

    for (auto k = 0; k < number_of_groups; ++k)
    {
        segments.push_back(std::round(normal_dist(generator)));
        chains.push_back(n0 * binomial_pmf(number_of_groups - 1, k));
    }

    if (!is_approx(n0, ranges::accumulate(chains, 0.0)))
    {
        std::cout << "My method failed\n";
        std::abort();
    }

    // Partition the range based on the average segments in order
    ranges::stable_partition(ranges::action::sort(segments),
                             [&](auto const N) { return N <= N_avg; });
}

std::vector<double> MicromechanicalElastomer::update_chains(
    std::vector<double> const& chains_old, double const time_step_size)
{
    using namespace ranges;

    return view::zip(segments, chains_old) | view::transform([&](auto const tpl) {
               auto const & [ N, n ] = tpl;
               return n / (1.0 + time_step_size * (1.0 - std::pow(1.0 - p_scission, N)));
           });
}

std::vector<double> MicromechanicalElastomer::compute_shear_moduli(
    std::vector<double> const& chains_new)
{
    return chains_new | ranges::view::transform([&](auto const& n) {
               return n * boltzmann_constant * temperature;
           });
}
}
