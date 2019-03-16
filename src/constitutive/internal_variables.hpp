
#pragma once

/// @file

#include "internal_variables_forward.hpp"

#include "numeric/dense_matrix.hpp"
#include "constitutive/variable_types.hpp"

#include <functional>
#include <map>
#include <vector>
#include <cstdint>

namespace neon
{
/// internal_variables stores the internal variables associated with the element
/// quadrature points.  These variables are duplicated and commited to memory
/// when the data is converged to avoid polluting the variable history in the
/// Newton-Raphson method.
template <typename SecondTensorType, typename FourthTensorType>
class internal_variables
{
public:
    /// A second order tensor type is a small matrix in tensor notation
    using second_tensor_type = SecondTensorType;
    /// A fourth order tensor type is a fixed size matrix in Voigt notation
    using fourth_tensor_type = FourthTensorType;

public:
    internal_variables(std::size_t const size) : m_size{size} {}

    /// Delete copy constructor to prevent references moving
    internal_variables(internal_variables const&) = delete;

    /// Delete assignment constructor to prevent references moving
    internal_variables& operator=(internal_variables const&) = delete;

    /// Implicitly defined move constructor
    internal_variables(internal_variables&&) = default;

    /// Add a number of scalar type variables to the object store
    template <typename... all_types>
    void add(variable::scalar const name, all_types const... names)
    {
        scalars[name].resize(m_size, 0.0);
        scalars_old[name].resize(m_size, 0.0);
        add(names...);
    }

    /// Add a number of tensor type variables to the object store
    template <typename... all_types>
    void add(variable::second const name, all_types... names)
    {
        second_order_tensors[name].resize(m_size, second_tensor_type::Zero());
        second_order_tensors_old[name].resize(m_size, second_tensor_type::Zero());
        add(names...);
    }

    /// Allocate scalars (defaulted to zeros)
    void add(variable::scalar const name, double const value = 0.0)
    {
        scalars[name].resize(m_size, value);
        scalars_old[name].resize(m_size, value);
    }

    /// Allocate second order tensors (defaulted to zeros)
    void add(variable::second const name)
    {
        second_order_tensors[name].resize(m_size, second_tensor_type::Zero());
        second_order_tensors_old[name].resize(m_size, second_tensor_type::Zero());
    }

    /// Allocate fourth order tensor (defaulted to zeros)
    void add(variable::fourth const name, fourth_tensor_type const m = fourth_tensor_type::Zero())
    {
        fourth_order_tensors[name].resize(m_size, m);
    }

    bool has(variable::scalar const name) const { return scalars.find(name) != scalars.end(); }

    bool has(variable::second const name) const
    {
        return second_order_tensors.find(name) != second_order_tensors.end();
    }

    bool has(variable::fourth const name) const
    {
        return fourth_order_tensors.find(name) != fourth_order_tensors.end();
    }

    /// Const access to the converged tensor variables
    std::vector<second_tensor_type> const& get_old(variable::second const name) const
    {
        return second_order_tensors_old.find(name)->second;
    }

    /// Const access to the converged scalar variables
    std::vector<double> const& get_old(variable::scalar const name) const
    {
        return scalars_old.find(name)->second;
    }

    /// Mutable access to the non-converged scalar variables
    std::vector<double>& get(variable::scalar const name)
    {
        if (!has(name))
        {
            throw std::domain_error("Scalar " + std::to_string(static_cast<int>(name))
                                    + " does not exist in the variable table");
        }
        return scalars.find(name)->second;
    }

    /// Mutable access to the non-converged second order tensor variables
    std::vector<second_tensor_type>& get(variable::second const name)
    {
        if (!has(name))
        {
            throw std::domain_error("Second order tensor " + std::to_string(static_cast<int>(name))
                                    + " does not exist in the variable table");
        }
        return second_order_tensors.find(name)->second;
    }

    /// Mutable access to the non-converged fourth order tensor variables
    std::vector<fourth_tensor_type>& get(variable::fourth const name)
    {
        if (!has(name))
        {
            throw std::domain_error("Fourth order tensor " + std::to_string(static_cast<int>(name))
                                    + " does not exist in the variable table");
        }
        return fourth_order_tensors.find(name)->second;
    }

    /// Mutable access to the non-converged scalar variables
    template <typename... scalar_types>
    auto get(variable::scalar const var0, scalar_types const... vars)
    {
        return std::make_tuple(std::ref(scalars.find(var0)->second),
                               std::ref(scalars.find(vars)->second)...);
    }

    /// Mutable access to the non-converged tensor variables
    template <typename... second_types>
    auto get(variable::second const var0, second_types const... vars)
    {
        return std::make_tuple(std::ref(second_order_tensors.find(var0)->second),
                               std::ref(second_order_tensors.find(vars)->second)...);
    }

    template <typename... fourth_types>
    auto get(variable::fourth const var0, fourth_types const... vars)
    {
        return std::make_tuple(std::ref(fourth_order_tensors.find(var0)->second),
                               std::ref(fourth_order_tensors.find(vars)->second)...);
    }

    /// Constant access to the non-converged scalar variables
    std::vector<double> const& get(variable::scalar const name) const
    {
        return scalars.find(name)->second;
    }

    /// Non-mutable access to the non-converged tensor variables
    std::vector<second_tensor_type> const& get(variable::second const name) const
    {
        return second_order_tensors.find(name)->second;
    }

    /// Non-mutable access to the non-converged matrix variables
    std::vector<fourth_tensor_type> const& get(variable::fourth const name) const
    {
        return fourth_order_tensors.find(name)->second;
    }

    /// Const access to the non-converged scalar variables
    template <typename... scalar_types>
    auto get(variable::scalar const var0, scalar_types const... vars) const
    {
        return std::make_tuple(std::cref(scalars.find(var0)->second),
                               std::cref(scalars.find(vars)->second)...);
    }

    /// Const access to the non-converged tensor variables
    template <typename... tensor_types>
    auto get(variable::second const var0, tensor_types const... vars) const
    {
        return std::make_tuple(std::cref(second_order_tensors.find(var0)->second),
                               std::cref(second_order_tensors.find(vars)->second)...);
    }

    template <typename... fourth_types>
    auto get(variable::fourth const var0, fourth_types const... vars) const
    {
        return std::make_tuple(std::cref(fourth_order_tensors.find(var0)->second),
                               std::cref(fourth_order_tensors.find(vars)->second)...);
    }

    /// Commit to history when iteration converges
    void commit()
    {
        scalars_old = scalars;
        second_order_tensors_old = second_order_tensors;
    }

    /// Revert to the old state when iteration doesn't converge
    void revert()
    {
        scalars = scalars_old;
        second_order_tensors = second_order_tensors_old;
    }

    /// \return Number of internal variables
    auto entries() const noexcept -> std::size_t { return m_size; }

protected:
    /// scalar history
    std::map<variable::scalar, std::vector<double>> scalars;
    /// old scalar history
    std::map<variable::scalar, std::vector<double>> scalars_old;

    /// second order tensors
    std::map<variable::second, std::vector<second_tensor_type>> second_order_tensors;
    /// old second order tensors
    std::map<variable::second, std::vector<second_tensor_type>> second_order_tensors_old;

    /// Fourth order tensors
    std::map<variable::fourth, std::vector<fourth_tensor_type>> fourth_order_tensors;

    std::size_t m_size;
};
}
