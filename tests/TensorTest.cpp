
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "numeric/Tensor.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Tensor decompositions")
{
    Matrix3 t = Matrix3::Random();

    auto const t_dev = deviatoric(t);
    auto const t_vol = volumetric(t);

    REQUIRE(((t_dev + t_vol) - t).norm() == Approx(0.0));
    REQUIRE(t_dev.trace() == Approx(0.0).margin(ZERO_MARGIN));
    REQUIRE(t_vol.trace() == Approx(t.trace()));
}
