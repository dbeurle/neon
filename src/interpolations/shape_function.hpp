
#pragma once

/// @file

#include <memory>
#include <utility>
#include <tuple>
#include <vector>

#include "numeric/dense_matrix.hpp"

namespace neon
{
/// shape_function is a base class for all polynomial based shape functions
template <typename... Spaces>
class shape_function
{
public:
    static constexpr auto spatial_dimension = sizeof...(Spaces);

    using coordinate_type = std::tuple<int, Spaces...>;

    /// Fix the size of the shape function derivative to the size of the quadrature points
    using value_type = std::pair<vector, matrixxd<sizeof...(Spaces)>>;

public:
    /// Construct the shape function by specifying
    /// \param node_count Number of nodes per element
    /// \param polynomial_order Highest polynomial order of shape function
    /// \param monomial_order Highest monomial order of shape function
    explicit shape_function(std::uint8_t const node_count,
                            std::uint8_t const polynomial_order,
                            std::uint8_t const monomial_order)
        : m_node_count{node_count}, m_polynomial_order{polynomial_order}, m_monomial_order{monomial_order}
    {
    }

    virtual ~shape_function() = default;

    /// \return Number of nodes in the element
    [[nodiscard]] auto number_of_nodes() const -> std::uint8_t { return m_node_count; }
    /// \return Highest polynomial order in interpolation function
    [[nodiscard]] auto polynomial_order() const -> std::uint8_t { return m_polynomial_order; }
    /// \return Highest monomial order in interpolation function
    [[nodiscard]] auto monomial_order() const -> std::uint8_t { return m_monomial_order; }

    /// Evaluate the shape functions at the natural coordinate
    [[nodiscard]] virtual auto evaluate(coordinate_type const&) const noexcept(false)
        -> value_type = 0;

    /// Evaluate the shape functions at multiple natural coordinates
    [[nodiscard]] auto evaluate(std::vector<coordinate_type> const& coordinates) const
        noexcept(false) -> std::vector<value_type>
    {
        // This trick works because the virtual dispatch will pickup the correct
        // evaluate() function call from above.  We are barred from instantiating
        // a shape_function object directly since we have a pure virtual method
        std::vector<value_type> values;
        values.reserve(coordinates.size());

        for (auto const& coordinate : coordinates)
        {
            values.emplace_back(this->evaluate(coordinate));
        }
        return values;
    }

    /// \return Natural element coordinates
    [[nodiscard]] auto local_coordinates() const noexcept -> std::vector<coordinate_type> const&
    {
        return m_local_coordinates;
    }

protected:
    /// Nodes per element
    std::uint8_t m_node_count{0};
    /// Highest order of polynomral term
    std::uint8_t m_polynomial_order{0};
    /// Highest order of monomial term
    std::uint8_t m_monomial_order{0};

    /// Natural element coordinates
    std::vector<coordinate_type> m_local_coordinates;
};

extern template class shape_function<double>;
extern template class shape_function<double, double>;
extern template class shape_function<double, double, double>;

using line_interpolation = shape_function<double>;
using surface_interpolation = shape_function<double, double>;
using volume_interpolation = shape_function<double, double, double>;
}
