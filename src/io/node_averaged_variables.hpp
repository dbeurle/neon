
#pragma once

#include "numeric/dense_matrix.hpp"

/// \file node_averaged_variable.hpp

namespace neon
{
/// Interpolate the internal variables to the mesh nodes and perform an
/// unweighted average.
template <typename mesh_type, typename enum_type>
vector average_internal_variable(mesh_type const& meshes,
                                 std::int32_t const variable_size,
                                 std::string const& variable_name,
                                 enum_type const variable_enum)
{
    vector nodal_averaged_value = vector::Zeros(mesh.geometry().size() * variable_size);
    vector insertion_count = vector::Zeros(mesh.geometry().size() * variable_size);

    // Add internal variables
    for (auto const& submesh : meshes)
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
}
