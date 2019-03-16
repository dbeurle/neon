
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"
#include "math/view.hpp"

#include <cstdint>
#include <string>
#include <memory>

namespace neon
{
/// Interpolate the internal variables to the mesh nodes and perform an
/// unweighted average.
template <typename MeshType, typename VariableType>
vector average_internal_variable(MeshType const& submeshes,
                                 std::int64_t const variable_size,
                                 std::string const& variable_name,
                                 VariableType const variable_enum)
{
    vector nodal_averaged_value = vector::Zero(variable_size);
    vector insertion_count = vector::Zero(variable_size);

    // Add internal variables
    for (auto const& submesh : submeshes)
    {
        if (!submesh.internal_variables().has(variable_enum))
        {
            throw std::domain_error("Internal variable " + variable_name + " does not exist in mesh");
        }

        auto const [value, count] = submesh.nodal_averaged_variable(variable_enum);

        nodal_averaged_value += value;

        insertion_count += count;
    }
    return nodal_averaged_value.cwiseQuotient(insertion_count);
}

template <class VariableType, class InternalVariableType, class IndicesType, std::size_t Size = 1>
auto nodal_averaged_variable(VariableType const name,
                             std::int32_t const component_count,
                             std::int64_t const node_count,
                             std::int64_t const element_count,
                             stride_view<> view,
                             IndicesType const& node_indices,
                             matrix const& extrapolation_matrix,
                             std::shared_ptr<InternalVariableType const> const& variables) noexcept(false)
    -> std::pair<vector, vector>
{
    vector count = vector::Zero(node_count * Size * Size);
    vector value = count;

    auto const& variable_list = variables->get(name);

    // vector format of values
    vector component(component_count);

    for (std::int64_t element{0}; element < element_count; ++element)
    {
        for (std::size_t i = 0; i < Size; ++i)
        {
            for (std::size_t j = 0; j < Size; ++j)
            {
                for (std::size_t l{0}; l < component_count; ++l)
                {
                    component(l) = variable_list[view(element, l)](i, j);
                }

                // Local extrapolation to the nodes
                vector const nodal_component = extrapolation_matrix * component;

                for (auto n = 0; n < nodal_component.rows(); n++)
                {
                    value(node_indices[n] * Size * Size + i * Size + j) += nodal_component(n);
                    count(node_indices[n] * Size * Size + i * Size + j) += 1.0;
                }
            }
        }
    }
    return {value, count};
}

}
