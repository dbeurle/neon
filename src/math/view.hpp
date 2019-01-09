
#pragma once

/// @file

#include <cstdint>
#include <type_traits>

namespace neon
{
/**
 * transform_expand_view takes a view of a contiguous input and output view
 * of a contiguous storage and applies an inline expansion according to the
 * ordering provided.
 *
 * An example of this transformation is given by
 *
 * Input: (3, 6, 7)
 * Order: {0, 1, 2}
 * Output: (9, 10, 11, 18, 19, 20, 21, 22, 23)
 *
 *
 * \param [in] input_view Contiguous view of input data to be transformed
 * \param [in out] output_view Contiguous view of input data to be transformed
 * \param order Storage e.g. {0, 1, 2}
 */
template <class InputView, class OutputView, typename OrderType>
inline void transform_expand_view(InputView input_view, OutputView output_view, OrderType const& order)
{
    using index_type = decltype(output_view.size());

    static_assert(std::is_integral<index_type>::value, "The size of the view needs to be an integer");

    static_assert(sizeof(typename InputView::value_type) >= sizeof(typename OutputView::value_type),
                  "The value type for the output range must be greater than or equal to input "
                  "range");

    for (index_type i{0}; i < input_view.size() * static_cast<index_type>(order.size()); ++i)
    {
        index_type const input_index = i / static_cast<index_type>(order.size());
        index_type const order_index = i % static_cast<index_type>(order.size());

        output_view[i] = input_view[input_index] * order.size() + order[order_index];
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
