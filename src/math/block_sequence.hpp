
#pragma once

#include <cstdint>

/// \file block_sequence.hpp

namespace neon
{
/// block_sequence provides the index mapping for Eigen views.  This functionality
/// is useful when we want to decode from the output of a linear solve where the
/// unknowns are encoded.
/// \tparam BlockSize Contiguous block of values to read before incrementing
/// \tparam IncrementSize The increment in the sequence
template <std::int64_t BlockSize = 1, std::int64_t IncrementSize = 1>
class block_sequence
{
public:
    using index_type = std::int64_t;

    static_assert(BlockSize <= IncrementSize,
                  "BlockSize must be less than or equal to increment_size");
    static_assert(IncrementSize > 0, "Increments must be greater than 0");
    static_assert(BlockSize > 0, "Blocks must be greater than 0");

public:
    /// Construct the sequence
    /// \param first The first element in the sequence
    /// \param size The total number of values in the sequence (<N)
    block_sequence(index_type const first, index_type const size) : first(first), m_size(size) {}

    /// \return Sequence size
    auto size() const noexcept { return m_size; }

    /// \return Given a sequence index provide the index into the original vector
    auto operator[](index_type const index) const noexcept
    {
        auto const quotient = index % BlockSize;
        auto const division = index / BlockSize;

        return first + IncrementSize * division + quotient;
    }

protected:
    index_type first, m_size;
};
}
