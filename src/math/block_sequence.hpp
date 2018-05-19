
#pragma once

#include <cstdint>

namespace neon
{
/// block_sequence provides the index mapping for Eigen views.  This functionality
/// is useful when we want to decode from the output of a linear solve where the
/// unknowns are encoded.
/// \tparam increment_size The increment in the sequence
/// \tparam block_size Contiguous values to read in
template <std::int64_t increment_size = 1, std::int64_t block_size = 1>
class block_sequence
{
public:
    using index_type = std::int64_t;

    static_assert(std::is_integral<index_type>::value, "Indices need to be integer types");
    static_assert(increment_size > 0, "Increments must be greater than 0");
    static_assert(block_size > 0, "Blocks must be greater than 0");

public:
    /// Construct the sequence
    /// \param first The first element in the sequence
    /// \param size The total number of values in the sequence (<N)
    block_sequence(index_type const first, index_type const size) : m_first(first), m_size(size) {}

    /// \return Sequence size
    auto size() const { return m_size; }

    /// \return Given a sequence index provide the index into the original vector
    auto operator[](index_type const index) const
    {
        auto const quotient = index % block_size;
        auto const remainder = index / block_size;

        return m_first + increment_size * remainder + quotient;
    }

protected:
    index_type m_first, m_size;
};
}
