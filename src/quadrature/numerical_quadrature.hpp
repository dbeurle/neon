
#pragma once

#include <tuple>
#include <vector>

namespace neon
{
/// numerical_quadrature is a variadic base class that defines a generic numerical
/// quadrature class that encapsulates the coordinates, weightings and a method
/// to perform the integration for a function that accepts a quadrature point
/// index.
template <typename... Spaces>
class numerical_quadrature
{
public:
    static constexpr auto spatial_dimension = sizeof...(Spaces);

    using coordinate_type = std::tuple<int, Spaces...>;

public:
    /// \return number of quadrature points
    [[nodiscard]] auto points() const noexcept { return m_weights.size(); }

    /// \return quadrature weights for this scheme
    [[nodiscard]] auto weights() const noexcept -> std::vector<double> const& { return m_weights; }

    /// \return list of tuples representing the index and the coordinates
    [[nodiscard]] auto coordinates() const noexcept -> auto const& { return m_coordinates; }

    /// \return Degree of polynomial exactly integrated
    [[nodiscard]] auto degree() const noexcept -> std::uint8_t { return m_degree; }

protected:
    /// Quadrature weightings
    std::vector<double> m_weights;
    /// Quadrature coordinates
    std::vector<coordinate_type> m_coordinates;

    std::uint8_t m_degree{0};
};

using line_quadrature = numerical_quadrature<double>;
using surface_quadrature = numerical_quadrature<double, double>;
using volume_quadrature = numerical_quadrature<double, double, double>;

extern template class numerical_quadrature<double>;
extern template class numerical_quadrature<double, double>;
extern template class numerical_quadrature<double, double, double>;

}
