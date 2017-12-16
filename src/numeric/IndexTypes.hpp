
#pragma once

#include <cstdint>
#include <vector>

namespace neon
{
using int32 = std::int32_t;
using int64 = std::int64_t;

using List = std::vector<int32>;

using local_indices = std::vector<std::int32_t>;
using nonlocal_indices = std::vector<std::int64_t>;
}
