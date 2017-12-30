
#pragma once

#include "numeric/IndexTypes.hpp"

#include <unordered_map>

namespace neon::boundary
{
class interprocess
{
public:
    interprocess(std::vector<int64> local_to_nonlocal,
                 std::unordered_map<int32, std::vector<int64>> process_interfaces)
        : local_to_nonlocal(local_to_nonlocal), process_interfaces(process_interfaces)
    {
    }

    auto shared_across() const { return process_interfaces.size(); }

    auto const& local_to_nonlocal_ordering() const { return local_to_nonlocal; }

    auto const& interfaces() const { return process_interfaces; }

protected:
    std::vector<int64> local_to_nonlocal;

    std::unordered_map<int32, std::vector<int64>> process_interfaces;
};
}
