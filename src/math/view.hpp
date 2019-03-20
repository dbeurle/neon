
#pragma once

/// @file

#include <cstdint>
#include <type_traits>

namespace neon
{
/**
 * transform_expand_view takes a view of a contiguous input and output iterator
 * of a contiguous storage and applies an inline expansion according to the
 * ordering provided.
 *
 * An example of this transformation is given by
 *
 * Input: (3, 6, 7)
 * Size: 3
 * Output: (9, 10, 11, 18, 19, 20, 21, 22, 23)
 *
 *
 * \param [in] input Contiguous iterator of input data to be transformed
 * \param [in] count A size to use
 * \param [in out] output Contiguous iterator of input data to be transformed (must be larger than
 * input)
 * \param [in] size to perform the inner loop to expand
 */
template <class InputIterator, class Size, class OutputIterator, class ExpandSize>
inline void transform_expand_n(InputIterator input,
                               Size const count,
                               OutputIterator output,
                               ExpandSize const size) noexcept
{
    static_assert(std::is_integral<Size>::value, "count must be an integer type");
    static_assert(std::is_integral<ExpandSize>::value, "size must be an integer type");

    for (Size outer_index = 0; outer_index < count; ++outer_index)
    {
        for (ExpandSize inner_index = 0; inner_index < size; ++inner_index)
        {
            *output++ = *input * size + inner_index;
        }
        ++input;
    }
}

/// A wrapper around the linear indexing into a compacted array using a uniform
/// stride index similar to a row-major ordering matrix in a one-dimensional
/// vector
template <class Index = std::size_t>
class stride_view
{
public:
    static_assert(std::is_integral<Index>::value, "Index must be an integer");

public:
    /// Construct with the stride size
    stride_view(Index const stride) : stride(stride) {}

    /// \return index into the linear storage
    auto operator()(Index const first_index, Index const second_index) const noexcept -> Index
    {
        return stride * first_index + second_index;
    }

private:
    Index stride;
};

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
