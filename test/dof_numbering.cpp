
#include <catch2/catch.hpp>

#include "numeric/index_types.hpp"
#include "math/view.hpp"

#include <memory>

using namespace neon;

TEST_CASE("index_view allocation")
{
    auto constexpr dof_index = 3;

    SECTION("One element")
    {
        indices nodal_connectivity(4, 1);
        nodal_connectivity << 0, 1, 2, 3;

        indices known_dof_list(4 * 3, 1);
        known_dof_list << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;

        indices computed_list(4 * 3, 1);
        transform_expand_n(nodal_connectivity(Eigen::all, 0).data(),
                           nodal_connectivity(Eigen::all, 0).size(),
                           computed_list(Eigen::all, 0).data(),
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
        transform_expand_n(nodal_connectivity(Eigen::all, 0).data(),
                           nodal_connectivity(Eigen::all, 0).size(),
                           computed_list(Eigen::all, 0).data(),
                           dof_index);

        transform_expand_n(nodal_connectivity(Eigen::all, 1).data(),
                           nodal_connectivity(Eigen::all, 1).size(),
                           computed_list(Eigen::all, 1).data(),
                           dof_index);

        REQUIRE((known_dof_list - computed_list).sum() == 0);
    }
}
