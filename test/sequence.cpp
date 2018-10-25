
#include <catch.hpp>

#include "math/view.hpp"

TEST_CASE("one sequence")
{
    neon::block_sequence<> sequence(0, 10);

    REQUIRE(sequence.size() == 10);

    REQUIRE(sequence[0] == 0);
    REQUIRE(sequence[1] == 1);
    REQUIRE(sequence[2] == 2);
    REQUIRE(sequence[10] == 10);
    REQUIRE(sequence[20] == 20);
    REQUIRE(sequence[30] == 30);
    REQUIRE(sequence[40] == 40);
    REQUIRE(sequence[50] == 50);
}
TEST_CASE("stride")
{
    {
        neon::stride_view<> view(1);

        REQUIRE(view(10, 0) == 10);
        REQUIRE(view(10, 2) == 12);
    }
    {
        neon::stride_view<> view(0);

        REQUIRE(view(10, 0) == 0);
        REQUIRE(view(10, 2) == 2);
    }
    {
        neon::stride_view<> view(2);

        REQUIRE(view(10, 0) == 20);
        REQUIRE(view(10, 2) == 22);
    }
}
