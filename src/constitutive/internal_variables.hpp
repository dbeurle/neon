
#pragma once

#include "internal_variables_forward.hpp"

#include "numeric/tensor_operations.hpp"

#include <functional>
#include <unordered_map>
#include <vector>
#include <cstdint>

namespace neon
{
class variable_view
{
public:
    variable_view() = default;

    variable_view(std::size_t const quadrature_points) : quadrature_points(quadrature_points) {}

    /** @return view index into the vector */
    std::size_t operator()(std::size_t const element, std::size_t const quadrature_point) const
    {
        return quadrature_points * element + quadrature_point;
    }

private:
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
    /** Spatial dimension (three, two or one dimension) */
    static auto constexpr r_n = rank2_dimension;
    /** Voigt dimension for the tensor conversion */
    static auto constexpr v_n = rank4_dimension;

    using scalar_type = double;

    /** A second order tensor type is a small matrix in tensor notation */
    using rank2tensor_type = Eigen::Matrix<scalar_type, rank2_dimension, rank2_dimension>;
    /** A fourth order tensor type is a fixed size matrix in Voigt notation */
    using rank4tensor_type = Eigen::Matrix<scalar_type, rank4_dimension, rank4_dimension>;
    static auto constexpr tensor_size = rank2_dimension * rank2_dimension;

public:
    enum class fourth : std::uint8_t { tangent_operator };

    /** Enumerations for internal variables that are tensor types */
    enum class second : std::uint8_t {
        /* Tensors for solid mechanics applications */
        // Stress measures
        Cauchy,
        Kirchhoff,
        PiolaKirchhoff1,
        PiolaKirchhoff2,
        // Deformation measures
        LinearisedStrain,
        LinearisedPlasticStrain,
        HenckyStrainElastic,
        DeformationGradient,
        DeformationGradientPlastic,
        DisplacementGradient,
        GreenLagrange,
        //  small_strain_J2_plasticity internal variables
        BackStress,
        KinematicHardening,
        /* Tensors for diffusion applications */
        Conductivity
    };

    enum class vector : std::uint8_t {
        Chains,
        Segments,
        ShearModuli,
        HeatFlux,
        deformation_history
    };

    enum class scalar : std::uint8_t {
        Chains,      // Chains for the micromechanical model
        Segments,    // Segments for the micromechanical model
        ShearModuli, // Shear moduli
        VonMisesStress,
        EffectivePlasticStrain,
        DetF0, // Reference Jacobian determinant
        DetF,  // Updated Jacobian determinant
        Damage,
        EnergyReleaseRate
    };

public:
    internal_variables(std::size_t const size) : size{size} {}

    /** Delete copy constructor to prevent references moving */
    internal_variables(internal_variables const&) = delete;

    /** Delete assignment constructor to prevent references moving */
    internal_variables& operator=(internal_variables const&) = delete;

    /** Implicitly defined move constructor */
    internal_variables(internal_variables&&) = default;

    /** Add a number of tensor type variables to the object store */
    template <typename... Variables>
    void add(second const name, Variables... names)
    {
        rank2tensors[name].resize(size, rank2tensor_type::Zero());
        rank2tensors_old[name].resize(size, rank2tensor_type::Zero());
        add(names...);
    }

    void add(second const name)
    {
        rank2tensors[name].resize(size, rank2tensor_type::Zero());
        rank2tensors_old[name].resize(size, rank2tensor_type::Zero());
    }

    /** Add a number of scalar type variables to the object store */
    template <typename... Variables>
    void add(scalar const name, Variables const... names)
    {
        scalars[name].resize(size, 0.0);
        scalars_old[name].resize(size, 0.0);
        add(names...);
    }

    void add(scalar const name)
    {
        scalars[name].resize(size, 0.0);
        scalars_old[name].resize(size, 0.0);
    }

    /** Allocate matrix internal variables with provided matrix */
    void add(internal_variables::fourth const name,
             rank4tensor_type const m = rank4tensor_type::Zero())
    {
        rank4tensors[name].resize(size, m);
    }

    bool has(scalar const name) const { return scalars.find(name) != scalars.end(); }

    bool has(second const name) const { return rank2tensors.find(name) != rank2tensors.end(); }

    bool has(fourth const name) const { return rank4tensors.find(name) != rank4tensors.end(); }

    /** Const access to the converged tensor variables */
    std::vector<rank2tensor_type> const& fetch_old(second const tensor_name) const
    {
        return rank2tensors_old.find(tensor_name)->second;
    }

    /** Const access to the converged scalar variables */
    std::vector<scalar_type> const& fetch_old(scalar const scalar_name) const
    {
        return scalars_old.find(scalar_name)->second;
    }

    /*-------------------------------------------------------------*
     *  Mutable access methods for unconverged internal variables  *
     *-------------------------------------------------------------*/

    /** Mutable access to the non-converged scalar variables */
    std::vector<scalar_type>& fetch(scalar const scalar_name)
    {
        return scalars.find(scalar_name)->second;
    }

    /** Mutable access to the non-converged tensor variables */
    std::vector<rank2tensor_type>& fetch(second tensor_name)
    {
        return rank2tensors.find(tensor_name)->second;
    }

    /** Mutable access to the non-converged matrix variables */
    std::vector<rank4tensor_type>& fetch(fourth const rank4_name)
    {
        return rank4tensors.find(rank4_name)->second;
    }

    /** Mutable access to the non-converged scalar variables */
    template <typename... scalar_types>
    auto fetch(scalar const var0, scalar const var1, scalar_types const... vars)
    {
        return std::make_tuple(std::ref(scalars.find(var0)->second),
                               std::ref(scalars.find(var1)->second),
                               std::ref(scalars.find(vars)->second)...);
    }

    /** Mutable access to the non-converged tensor variables */
    template <typename... TensorTps>
    auto fetch(second var0, second var1, TensorTps... vars)
    {
        return std::make_tuple(std::ref(rank2tensors.find(var0)->second),
                               std::ref(rank2tensors.find(var1)->second),
                               std::ref(rank2tensors.find(vars)->second)...);
    }

    template <typename... rank4_types>
    auto fetch(fourth const var0, fourth const var1, rank4_types const... vars)
    {
        return std::make_tuple(std::ref(rank4tensors.find(var0)->second),
                               std::ref(rank4tensors.find(var1)->second),
                               std::ref(rank4tensors.find(vars)->second)...);
    }

    /*-------------------------------------------------------------*
     * Constant access methods for unconverged internal variables  *
     *-------------------------------------------------------------*/

    /** Constant access to the non-converged scalar variables */
    std::vector<scalar_type> const& fetch(scalar const scalar_name) const
    {
        return scalars.find(scalar_name)->second;
    }

    /** Non-mutable access to the non-converged tensor variables */
    std::vector<rank2tensor_type> const& fetch(second const tensor_name) const
    {
        return rank2tensors.find(tensor_name)->second;
    }

    /** Non-mutable access to the non-converged matrix variables */
    std::vector<rank4tensor_type> const& fetch(fourth const rank4_name) const
    {
        return rank4tensors.find(rank4_name)->second;
    }

    /** Const access to the non-converged scalar variables */
    template <typename... ScalarTps>
    auto fetch(scalar const var0, scalar const var1, ScalarTps const... vars) const
    {
        return std::make_tuple(std::cref(scalars.find(var0)->second),
                               std::cref(scalars.find(var1)->second),
                               std::cref(scalars.find(vars)->second)...);
    }

    /** Const access to the non-converged tensor variables */
    template <typename... TensorTps>
    auto fetch(second var0, second var1, TensorTps... vars) const
    {
        return std::make_tuple(std::cref(rank2tensors.find(var0)->second),
                               std::cref(rank2tensors.find(var1)->second),
                               std::cref(rank2tensors.find(vars)->second)...);
    }

    template <typename... rank4_types>
    auto fetch(fourth const var0, fourth const var1, rank4_types... vars) const
    {
        return std::make_tuple(std::cref(rank4tensors.find(var0)->second),
                               std::cref(rank4tensors.find(var1)->second),
                               std::cref(rank4tensors.find(vars)->second)...);
    }

    /** Commit to history when iteration converges */
    void commit()
    {
        rank2tensors_old = rank2tensors;
        scalars_old = scalars;
    }

    /** Revert to the old state when iteration doesn't converge */
    void revert()
    {
        rank2tensors = rank2tensors_old;
        scalars = scalars_old;
    }

    auto entries() const { return size; }

protected:
    // These state variables are committed and reverted depending on the outer
    // simulation loop.  If a nonlinear iteration does not converge then revert
    // the state back to the previous state.  The rank2tensors and scalars fields
    // are the 'unstable' variables and the *Old are the stable variants
    std::unordered_map<second, std::vector<rank2tensor_type>> rank2tensors, rank2tensors_old;

    std::unordered_map<scalar, std::vector<scalar_type>> scalars, scalars_old;

    std::unordered_map<fourth, std::vector<rank4tensor_type>> rank4tensors;

    std::size_t size;
};
}
