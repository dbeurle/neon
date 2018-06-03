
#pragma once

#include "internal_variables_forward.hpp"

#include "numeric/dense_matrix.hpp"

#include <functional>
#include <unordered_map>
#include <vector>
#include <cstdint>

namespace neon
{
/// variable_view is a tiny wrapper around the linear indexing into the
/// internal variable storage given an element quadrature point and an element
/// number.  \sa internal_variables
class variable_view
{
public:
    variable_view() = default;

    /// Construct with number of quadrature points per element
    variable_view(std::size_t const quadrature_points) : quadrature_points(quadrature_points) {}

    /// \return view index into the vector
    std::size_t operator()(std::size_t const element, std::size_t const quadrature_point) const
        noexcept
    {
        return quadrature_points * element + quadrature_point;
    }

private:
    /// Quadrature points per element
    std::size_t quadrature_points{0};
};

/// internal_variables stores the internal variables associated with the element
/// quadrature points.  These variables are duplicated and commited to memory
/// when the data is converged to avoid polluting the variable history in the
/// Newton-Raphson method.
template <int rank2_dimension, int rank4_dimension>
class internal_variables
{
public:
    /// Spatial dimension (three, two or one dimension)
    static auto constexpr r_n = rank2_dimension;
    /// Voigt dimension for the tensor conversion
    static auto constexpr v_n = rank4_dimension;

    using scalar_type = double;

    /// A second order tensor type is a small matrix in tensor notation
    using second_tensor_type = Eigen::Matrix<scalar_type, rank2_dimension, rank2_dimension>;
    /// A fourth order tensor type is a fixed size matrix in Voigt notation
    using fourth_tensor_type = Eigen::Matrix<scalar_type, rank4_dimension, rank4_dimension>;

    static auto constexpr tensor_size = rank2_dimension * rank2_dimension;

public:
    /// Fourth order tensor types
    enum class fourth : std::uint8_t {
        /// Material tangent operator
        tangent_operator
    };

    /// Second order tensor internal variables types
    enum class second : std::uint8_t {
        /// Cauchy stress
        cauchy_stress,
        /// Kirchhoff stress
        Kirchhoff,
        /// First Piola-Kirchhoff stress (PK1)
        PiolaKirchhoff1,
        /// Second Piola-Kirchhoff stress (PK2)
        PiolaKirchhoff2,
        /// Linearised (infinitesimal) total strain
        LinearisedStrain,
        /// Linearised (infinitesimal) plastic strain
        LinearisedPlasticStrain,
        /// Hencky elastic strain
        HenckyStrainElastic,
        /// Deformation gradient (F)
        DeformationGradient,
        /// Plastic deformation gradient (Fp)
        DeformationGradientPlastic,
        /// Displacement gradient (H)
        DisplacementGradient,
        /// Green-Lagrange strain (E)
        GreenLagrange,
        /// Back stress (for hardening)
        BackStress,
        /// Kinematic hardening
        KinematicHardening,
        /// Conductivity tensor
        Conductivity,
        /// Beam bending stiffness
        bending_stiffness,
        /// Beam shear stiffness
        shear_stiffness,
        /// Secondary network accumulated integral
        accumulated_ageing_integral,
        /// Last evaluated integrand for accumulated integral
        ageing_integrand
    };

    /// Names for vector values
    enum class vector : std::uint8_t { HeatFlux, deformation_history };

    /// Names for scalar values
    enum class scalar : std::uint8_t {
        /// Active chains per unit volume
        active_shear_modulus,
        /// Inactive chains per unit volume
        inactive_shear_modulus,
        /// Active number of segments per chain
        active_segments,
        /// Inactive number of segments per chain
        inactive_segments,
        /// Reduction factor for network decay
        reduction_factor,
        /// Von Mises equivalent stress
        VonMisesStress,
        /// Accumulated (effective) plastic strain
        EffectivePlasticStrain,
        /// Reference Jacobian determinant
        DetF0,
        /// Updated Jacobian determinant
        DetF,
        /// Scalar damage variable
        Damage,
        /// Scalar energy release rate
        EnergyReleaseRate,
        /// Beam torsional stiffness
        torsional_stiffness,
        /// Beam axial stiffness
        axial_stiffness,
        /// Beam section shear area (one)
        shear_area_1,
        /// Beam section shear area (two)
        shear_area_2,
        /// Section cross sectional area
        cross_sectional_area,
        /// Beam second moment of area (one)
        second_moment_area_1,
        /// Beam second moment of area (two)
        second_moment_area_2
    };

public:
    internal_variables(std::size_t const size) : size{size} {}

    /// Delete copy constructor to prevent references moving
    internal_variables(internal_variables const&) = delete;

    /// Delete assignment constructor to prevent references moving
    internal_variables& operator=(internal_variables const&) = delete;

    /// Implicitly defined move constructor
    internal_variables(internal_variables&&) = default;

    /// Add a number of tensor type variables to the object store
    template <typename... all_types>
    void add(second const name, all_types... names)
    {
        second_order_tensors[name].resize(size, second_tensor_type::Zero());
        second_order_tensors_old[name].resize(size, second_tensor_type::Zero());
        add(names...);
    }

    /// Add a number of scalar type variables to the object store
    template <typename... all_types>
    void add(scalar const name, all_types const... names)
    {
        scalars[name].resize(size, 0.0);
        scalars_old[name].resize(size, 0.0);
        add(names...);
    }

    /// Allocate scalars (defaulted to zeros)
    void add(scalar const name, double const value = 0.0)
    {
        scalars[name].resize(size, value);
        scalars_old[name].resize(size, value);
    }

    /// Allocate second order tensors (defaulted to zeros)
    void add(second const name)
    {
        second_order_tensors[name].resize(size, second_tensor_type::Zero());
        second_order_tensors_old[name].resize(size, second_tensor_type::Zero());
    }

    /// Allocate fourth order tensor (defaulted to zeros)
    void add(internal_variables::fourth const name,
             fourth_tensor_type const m = fourth_tensor_type::Zero())
    {
        fourth_order_tensors[name].resize(size, m);
    }

    bool has(scalar const name) const { return scalars.find(name) != scalars.end(); }

    bool has(second const name) const
    {
        return second_order_tensors.find(name) != second_order_tensors.end();
    }

    bool has(fourth const name) const
    {
        return fourth_order_tensors.find(name) != fourth_order_tensors.end();
    }

    /// Const access to the converged tensor variables
    std::vector<second_tensor_type> const& get_old(second const name) const
    {
        return second_order_tensors_old.find(name)->second;
    }

    /// Const access to the converged scalar variables
    std::vector<scalar_type> const& get_old(scalar const name) const
    {
        return scalars_old.find(name)->second;
    }

    /// Mutable access to the non-converged scalar variables
    std::vector<scalar_type>& get(scalar const name)
    {
        if (!has(name))
        {
            throw std::domain_error("Scalar " + std::to_string(static_cast<int>(name))
                                    + " does not exist in the std::unordered_map");
        }
        return scalars.find(name)->second;
    }

    /// Mutable access to the non-converged second order tensor variables
    std::vector<second_tensor_type>& get(second const name)
    {
        if (!has(name))
        {
            throw std::domain_error("Second order tensor " + std::to_string(static_cast<int>(name))
                                    + " does not exist in the std::unordered_map");
        }
        return second_order_tensors.find(name)->second;
    }

    /// Mutable access to the non-converged fourth order tensor variables
    std::vector<fourth_tensor_type>& get(fourth const name)
    {
        if (!has(name))
        {
            throw std::domain_error("Fourth order tensor " + std::to_string(static_cast<int>(name))
                                    + " does not exist in the std::unordered_map");
        }
        return fourth_order_tensors.find(name)->second;
    }

    /// Mutable access to the non-converged scalar variables
    template <typename... scalar_types>
    auto get(scalar const var0, scalar_types const... vars)
    {
        return std::make_tuple(std::ref(scalars.find(var0)->second),
                               std::ref(scalars.find(vars)->second)...);
    }

    /// Mutable access to the non-converged tensor variables
    template <typename... second_types>
    auto get(second const var0, second_types const... vars)
    {
        return std::make_tuple(std::ref(second_order_tensors.find(var0)->second),
                               std::ref(second_order_tensors.find(vars)->second)...);
    }

    template <typename... fourth_types>
    auto get(fourth const var0, fourth_types const... vars)
    {
        return std::make_tuple(std::ref(fourth_order_tensors.find(var0)->second),
                               std::ref(fourth_order_tensors.find(vars)->second)...);
    }

    /// Constant access to the non-converged scalar variables
    std::vector<scalar_type> const& get(scalar const name) const
    {
        return scalars.find(name)->second;
    }

    /// Non-mutable access to the non-converged tensor variables
    std::vector<second_tensor_type> const& get(second const name) const
    {
        return second_order_tensors.find(name)->second;
    }

    /// Non-mutable access to the non-converged matrix variables
    std::vector<fourth_tensor_type> const& get(fourth const name) const
    {
        return fourth_order_tensors.find(name)->second;
    }

    /// Const access to the non-converged scalar variables
    template <typename... scalar_types>
    auto get(scalar const var0, scalar_types const... vars) const
    {
        return std::make_tuple(std::cref(scalars.find(var0)->second),
                               std::cref(scalars.find(vars)->second)...);
    }

    /// Const access to the non-converged tensor variables
    template <typename... tensor_types>
    auto get(second const var0, tensor_types const... vars) const
    {
        return std::make_tuple(std::cref(second_order_tensors.find(var0)->second),
                               std::cref(second_order_tensors.find(vars)->second)...);
    }

    template <typename... fourth_types>
    auto get(fourth const var0, fourth_types const... vars) const
    {
        return std::make_tuple(std::cref(fourth_order_tensors.find(var0)->second),
                               std::cref(fourth_order_tensors.find(vars)->second)...);
    }

    /// Commit to history when iteration converges
    void commit()
    {
        second_order_tensors_old = second_order_tensors;
        scalars_old = scalars;
    }

    /// Revert to the old state when iteration doesn't converge
    void revert()
    {
        second_order_tensors = second_order_tensors_old;
        scalars = scalars_old;
    }

    /// \return Number of internal variables
    auto entries() const noexcept { return size; }

protected:
    // These state variables are committed and reverted depending on the outer
    // simulation loop.  If a nonlinear iteration does not converge then revert
    // the state back to the previous state.  The second_order_tensors and scalars fields
    // are the 'unstable' variables and the *old are the stable variants
    std::unordered_map<second, std::vector<second_tensor_type>> second_order_tensors,
        second_order_tensors_old;

    std::unordered_map<scalar, std::vector<scalar_type>> scalars, scalars_old;

    std::unordered_map<fourth, std::vector<fourth_tensor_type>> fourth_order_tensors;

    std::size_t size;
};
}
