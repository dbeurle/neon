
#include <catch.hpp>

#include "numeric/index_types.hpp"
#include "math/transform_expand.hpp"

#include <range/v3/view/set_algorithm.hpp>
#include <memory>

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

std::array<std::int32_t, 3> dof_index = {0, 1, 2};

TEST_CASE("index_view allocation")
{
    SECTION("One element")
    {
        indices nodal_connectivity(4, 1);
        nodal_connectivity << 0, 1, 2, 3;

        indices known_dof_list(4 * 3, 1);
        known_dof_list << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;

        indices computed_list(4 * 3, 1);
        transform_expand_view(nodal_connectivity(Eigen::placeholders::all, 0),
                              computed_list(Eigen::placeholders::all, 0),
                              dof_index);

        REQUIRE((known_dof_list - computed_list).sum() == 0);
    }
    SECTION("Two elements")
    {
        indices nodal_connectivity(2, 4);
        nodal_connectivity << 0, 1, 2, 3, 4, 2, 1, 5;
        nodal_connectivity.transposeInPlace();

        indices known_dof_list(2, 4 * 3);
        known_dof_list << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, //
            12, 13, 14, 6, 7, 8, 3, 4, 5, 15, 16, 17;
        known_dof_list.transposeInPlace();

        indices computed_list(4 * 3, 2);
        transform_expand_view(nodal_connectivity(Eigen::placeholders::all, 0),
                              computed_list(Eigen::placeholders::all, 0),
                              dof_index);

        transform_expand_view(nodal_connectivity(Eigen::placeholders::all, 1),
                              computed_list(Eigen::placeholders::all, 1),
                              dof_index);

        REQUIRE((known_dof_list - computed_list).sum() == 0);
    }
}
