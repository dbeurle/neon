
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"
#include "numeric/dense_matrix.hpp"

namespace neon
{
/// unit_sphere_quadrature implements the symmetric quadrature rules
class unit_sphere_quadrature : public volume_quadrature
{
public:
    /// Quadrature nomenclature from \cite Ehret2010
    enum class scheme { BO21, BO33, BO61, FM900 };

public:
    /// Construct the unit_sphere_quadrature using the available scheme
    unit_sphere_quadrature(scheme const rule);

    template <typename MatrixType, typename Callable>
    MatrixType integrate(MatrixType integrand, Callable&& f) const
    {
        auto const indices = m_weights.size();
        for (std::size_t index{0}; index < indices; ++index)
        {
            integrand += f(m_values[index], index) * m_weights[index];
        }
        return integrand;
    }

    template <typename Callable>
    double integrate(double integrand, Callable&& f) const
    {
        auto const indices = m_weights.size();
        for (std::size_t index{0}; index < indices; ++index)
        {
            integrand += f(m_values[index], index) * m_weights[index];
        }
        return integrand;
    }

    template <typename Callable>
    void for_each(Callable&& f) const
    {
        auto const indices = m_weights.size();
        for (std::size_t index{0}; index < indices; ++index)
        {
            f(m_values[index], index);
        }
    }

private:
    std::vector<std::pair<vector3, matrix3>> m_values;
};
}
