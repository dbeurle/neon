
#pragma once

#include <cstdint>
#include <utility>
#include <Eigen/Core>

namespace neon
{
using indices = Eigen::Array<std::int32_t, Eigen::Dynamic, Eigen::Dynamic>;

/// Type alias for whatever type is returned from these views
using index_view = decltype(std::declval<const indices>()(Eigen::placeholders::all, 0l));
}
