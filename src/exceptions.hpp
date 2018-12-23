
#pragma once

/// @file

#include <stdexcept>

namespace neon
{
/// computational_error should be thrown when something is wrong with the
/// computation that can potentially be corrected with a change in load step.
/// These occur during the linear solution phase, with element distortion and
/// violations of internal variable updates as examples
struct computational_error : public std::domain_error
{
    using std::domain_error::domain_error;
};

/// Exception type for computation error in CUDA kernel calls
struct cuda_error : public std::domain_error
{
    using std::domain_error::domain_error;
};
}
