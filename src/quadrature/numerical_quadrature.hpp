
#pragma once

#include <tuple>
#include <vector>

#include "numeric/dense_matrix.hpp"

namespace neon
{
/// numerical_quadrature is a variadic class that defines a generic numerical
/// quadrature class that encapsulates the coordinates, weightings and a method
/// to perform the integration for a function that accepts a quadrature point
/// index.
template <typename... Spaces>
class numerical_quadrature
{
public:
    static constexpr std::size_t spatial_dimension = sizeof...(Spaces);

    using coordinate_type = std::tuple<int, Spaces...>;

    /// Fix the size of the shape function derivative to the size of the quadrature points
    using value_type = std::tuple<vector, matrixxd<spatial_dimension>>;

public:
    /// Perform the numerical integration of a lambda function.
    /// \param integral - Initial value for the numerical integration
    /// \param f - A lambda function that accepts an femValue and quadrature point
    /// \return The numerically integrated matrix
    template <typename MatrixType, typename Callable>
    [[nodiscard]] MatrixType integrate(MatrixType operand, Callable&& f) const noexcept
    {
        for (std::size_t l{0}; l < points(); ++l)
        {
            operand.noalias() += f(values[l], l) * m_weights[l];
        }
        return operand;
    }

    /// Perform the numerical integration of a lambda function.
    /// \param integral - Initial value for the numerical integration
    /// \param f - A lambda function that accepts an femValue and quadrature point
    /// \return The numerically integrated scalar
    template <typename Callable>
    [[nodiscard]] double integrate(double integral, Callable&& f) const noexcept
    {
        for (std::size_t l{0}; l < points(); ++l)
        {
            integral += f(values[l], l) * m_weights[l];
        }
        return integral;
    }

    /// Perform the numerical integration of a lambda function.
    /// \param integral Value for the numerical integration (accumulated into)
    /// \param f A lambda function that accepts an femValue and quadrature point
    /// \return The numerically integrated matrix
    template <typename Callable>
    void integrate_inplace(matrix& integral, Callable&& f) const noexcept
    {
        for (std::size_t l{0}; l < points(); ++l)
        {
            integral.noalias() += f(values[l], l) * m_weights[l];
        }
    }

    /// Perform the numerical integration of a lambda function.
    /// \param integral Value for the numerical integration (accumulated into)
    /// \param f A lambda function that accepts an femValue and quadrature point
    /// \return The numerically integrated matrix
    template <typename MatrixType, typename Callable>
    void integrate_inplace(Eigen::Map<MatrixType> integral, Callable&& f) const noexcept
    {
        for (std::size_t index{0}; index < points(); ++index)
        {
            integral.noalias() += f(values[index], index) * m_weights[index];
        }
    }

    /// Evaluate the function \p function at each integration point
    /// \tparam A callable type
    template <typename Callable>
    void for_each(Callable&& function) const
    {
        for (std::size_t l{0}; l < points(); ++l)
        {
            function(values[l], l);
        }
    }

    /// Evaluate a shape function and derivatives for the quadrature points
    /// Throws std::out_of_memory
    /// \param f - A lambda function that accepts a quadrature coordinate tuple
    template <typename Callable>
    void evaluate(Callable&& f)
    {
        values.clear();
        values.reserve(points());
        for (auto const& coordinate : m_coordinates)
        {
            values.emplace_back(f(coordinate));
        }
    }

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

    /// Shape functions and their derivatives evaluated at the quadrature points
    std::vector<value_type> values;

    std::uint8_t m_degree{0};
};

extern template class numerical_quadrature<double>;
extern template class numerical_quadrature<double, double>;
extern template class numerical_quadrature<double, double, double>;

using line_quadrature = numerical_quadrature<double>;
using surface_quadrature = numerical_quadrature<double, double>;
using volume_quadrature = numerical_quadrature<double, double, double>;
}
