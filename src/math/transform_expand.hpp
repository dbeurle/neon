
#include <type_traits>

namespace neon
{
/**
 * transform_expand_view_n takes a view of a contiguous input and output view
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
template <class input_view_type, class output_view_type, typename order_type>
void transform_expand_view(input_view_type input_view,
                           output_view_type output_view,
                           order_type const& order)
{
    using index_type = decltype(output_view.size());

    static_assert(std::is_integral<index_type>::value, "The size of the view needs to be an integer");

    static_assert(sizeof(typename input_view_type::value_type)
                      >= sizeof(typename output_view_type::value_type),
                  "The value type for the output range must be greater than or equal to input "
                  "range");

    for (index_type i{0}; i < input_view.size() * order.size(); ++i)
    {
        index_type const input_index = i / static_cast<index_type>(order.size());
        index_type const order_index = i % static_cast<index_type>(order.size());

        output_view[i] = input_view[input_index] * order.size() + order[order_index];
    }
}
}
