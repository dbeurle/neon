
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/DofAllocator.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include <memory>

using namespace neon;
using namespace ranges;

TEST_CASE("Dof List Allocation", "[DofAllocator]")
{
    SECTION("One element")
    {
        std::vector<List> const nodal_connectivity = {{0, 1, 2, 3}};

        std::vector<List> const known_dof_list{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};

        std::vector<List> const computed_list = allocate_dof_list(3, nodal_connectivity);

        REQUIRE(view::set_difference(computed_list.at(0), known_dof_list.at(0)).empty());
    }
    SECTION("Two elements")
    {
        std::vector<List> const nodal_connectivity = {{0, 1, 2, 3}, {4, 2, 1, 5}};

        std::vector<List> const known_dof_list{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                                               {12, 13, 14, 6, 7, 8, 3, 4, 5, 15, 16, 17}};

        std::vector<List> const computed_list = allocate_dof_list(3, nodal_connectivity);

        REQUIRE(view::set_difference(computed_list.at(0), known_dof_list.at(0)).empty());
        REQUIRE(view::set_difference(computed_list.at(1), known_dof_list.at(1)).empty());
    }
}
