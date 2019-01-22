
#include "mesh/projection/recovery.hpp"

namespace neon
{
local_extrapolation::local_extrapolation(matrix const& shape_function,
                                         matrix const& local_nodal_coordinates,
                                         matrix const& local_quadrature_coordinates)
{
    this->compute_extrapolation_matrix(shape_function,
                                       local_nodal_coordinates,
                                       local_quadrature_coordinates);
}

auto local_extrapolation::project(std::vector<double> const& variables,
                                  stride_view<> const& stride,
                                  indices const& node_indices,
                                  std::int64_t const node_count) const noexcept -> value_type
{
    vector projections = vector::Zero(node_count);
    std::vector<std::uint16_t> counts(node_count, 0);

    // vector format of values
    vector component = vector::Zero(stride.size());

    auto const elements = node_indices.cols();

    for (std::int64_t element{0}; element < elements; ++element)
    {
        // Assemble these into the global value vector
        auto const& node_list = node_indices(element, Eigen::all);

        for (std::size_t index{0}; index < stride.size(); ++index)
        {
            component(index) = variables[stride(element, index)];
        }

        projections(node_list) += m_extrapolation * component;

        for (std::size_t index = 0; index < stride.size(); index++)
        {
            ++counts[node_list(index)];
        }
    }
    return {projections, counts};
}

void local_extrapolation::compute_extrapolation_matrix(matrix const& shape_function,
                                                       matrix const& local_nodal_coordinates,
                                                       matrix const& local_quadrature_coordinates)
{
    // Take short names for consistency with the algorithm presented in the paper

    // Narrowing conversion but rows is expected to be greater than zero
    std::size_t const n = local_nodal_coordinates.rows();
    std::size_t const m = local_quadrature_coordinates.rows();

    auto const& xi = local_nodal_coordinates;
    auto const& xi_hat = local_quadrature_coordinates;

    if (m == n)
    {
        m_extrapolation = shape_function.inverse();
        return;
    }

    if (m > n)
    {
        m_extrapolation = (shape_function.transpose() * shape_function).inverse()
                          * shape_function.transpose();
        return;
    }

    matrix const N_plus = shape_function.transpose()
                          * (shape_function * shape_function.transpose()).inverse();

    // Number of quadrature points are less than the number of nodes
    auto const xi_hat_plus = xi_hat.transpose() * (xi_hat * xi_hat.transpose()).inverse();

    m_extrapolation = N_plus * (matrix::Identity(m, m) - xi_hat * xi_hat_plus) + xi * xi_hat_plus;
}
}
