
#include <catch.hpp>

#include "math/block_sequence.hpp"

TEST_CASE("one sequence")
{
    neon::block_sequence<> sequence(0, 10);
}
