
#pragma once

#include "InternalVariablesForwards.hpp"

#include "numeric/tensor_operations.hpp"

#include <functional>
#include <unordered_map>
#include <vector>

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

/**
 * InternalVariables stores the internal variables associated with the element
 * quadrature points.  These variables are duplicated and commited to memory
 * when the data is converged to avoid polluting the variable history in the
 * Newton-Raphson method.
 */
template <int rank2_dimension, int rank4_dimension>
class InternalVariables
{
public:
    /** Spatial dimension (three, two or one dimension) */
    static auto constexpr r_n = rank2_dimension;

    /** Voigt dimension for the tensor conversion */
    static auto constexpr v_n = rank4_dimension;

    // Type aliases
    using scalar_type = double;

    /** A second order tensor type is a small matrix in tensor notation */
    using rank2tensor_type = Eigen::Matrix<scalar_type, rank2_dimension, rank2_dimension>;

    /** A fourth order tensor type is a fixed size matrix in Voigt notation */
    using rank4tensor_type = Eigen::Matrix<scalar_type, rank4_dimension, rank4_dimension>;

    static auto constexpr tensor_size = rank2_dimension * rank2_dimension;

public:
    enum class rank4 { tangent_operator };

    /** Enumerations for internal variables that are tensor types */
    enum class Tensor {
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

    enum class vector { Chains, Segments, ShearModuli, HeatFlux };

    enum class Scalar {
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
    InternalVariables(std::size_t const size) : size(size) {}

    /** Delete copy constructor to prevent references moving */
    InternalVariables(InternalVariables const&) = delete;

    /** Delete assignment constructor to prevent references moving */
    InternalVariables& operator=(InternalVariables const&) = delete;

    /** Implicitly defined move constructor */
    InternalVariables(InternalVariables&&) = default;

    /** Add a number of tensor type variables to the object store */
    template <typename... Variables>
    void add(Tensor const name, Variables... names)
    {
        rank2tensors[name].resize(size, rank2tensor_type::Zero());
        rank2tensors_old[name].resize(size, rank2tensor_type::Zero());
        add(names...);
    }

    void add(Tensor const name)
    {
        rank2tensors[name].resize(size, rank2tensor_type::Zero());
        rank2tensors_old[name].resize(size, rank2tensor_type::Zero());
    }

    /** Add a number of scalar type variables to the object store */
    template <typename... Variables>
    void add(Scalar const name, Variables... names)
    {
        scalars[name].resize(size, 0.0);
        scalars_old[name].resize(size, 0.0);
        add(names...);
    }

    void add(Scalar const name)
    {
        scalars[name].resize(size, 0.0);
        scalars_old[name].resize(size, 0.0);
    }

    /** Allocate matrix internal variables with provided matrix */
    void add(InternalVariables::rank4 const name, rank4tensor_type const m = rank4tensor_type::Zero())
    {
        rank4tensors[name].resize(size, m);
    }

    bool has(Scalar const name) const { return scalars.find(name) != scalars.end(); }
    bool has(Tensor const name) const { return rank2tensors.find(name) != rank2tensors.end(); }
    bool has(rank4 const name) const { return rank4tensors.find(name) != rank4tensors.end(); }

    /** Const access to the converged tensor variables */
    std::vector<rank2tensor_type> const& fetch_old(Tensor const tensorType) const
    {
        return rank2tensors_old.find(tensorType)->second;
    }

    /** Const access to the converged scalar variables */
    std::vector<scalar_type> const& fetch_old(Scalar const scalarType) const
    {
        return scalars_old.find(scalarType)->second;
    }

    /*-------------------------------------------------------------*
     *  Mutable access methods for unconverged internal variables  *
     *-------------------------------------------------------------*/

    /** Mutable access to the non-converged scalar variables */
    std::vector<scalar_type>& fetch(Scalar scalarType) { return scalars.find(scalarType)->second; }

    /** Mutable access to the non-converged tensor variables */
    std::vector<rank2tensor_type>& fetch(Tensor tensorType)
    {
        return rank2tensors.find(tensorType)->second;
    }

    /** Mutable access to the non-converged matrix variables */
    std::vector<rank4tensor_type>& fetch(rank4 const matrixType)
    {
        return rank4tensors.find(matrixType)->second;
    }

    /** Mutable access to the non-converged scalar variables */
    template <typename... ScalarTps>
    auto fetch(Scalar var0, Scalar var1, ScalarTps... vars)
    {
        return std::make_tuple(std::ref(scalars.find(var0)->second),
                               std::ref(scalars.find(var1)->second),
                               std::ref(scalars.find(vars)->second)...);
    }

    /** Mutable access to the non-converged tensor variables */
    template <typename... TensorTps>
    auto fetch(Tensor var0, Tensor var1, TensorTps... vars)
    {
        return std::make_tuple(std::ref(rank2tensors.find(var0)->second),
                               std::ref(rank2tensors.find(var1)->second),
                               std::ref(rank2tensors.find(vars)->second)...);
    }

    template <typename... rank4_types>
    auto fetch(rank4 var0, rank4 var1, rank4_types... vars)
    {
        return std::make_tuple(std::ref(rank4tensors.find(var0)->second),
                               std::ref(rank4tensors.find(var1)->second),
                               std::ref(rank4tensors.find(vars)->second)...);
    }

    /*-------------------------------------------------------------*
     * Constant access methods for unconverged internal variables  *
     *-------------------------------------------------------------*/

    /** Constant access to the non-converged scalar variables */
    std::vector<scalar_type> const& fetch(Scalar const scalarType) const
    {
        return scalars.find(scalarType)->second;
    }

    /** Non-mutable access to the non-converged tensor variables */
    std::vector<rank2tensor_type> const& fetch(Tensor const tensorType) const
    {
        return rank2tensors.find(tensorType)->second;
    }

    /** Non-mutable access to the non-converged matrix variables */
    std::vector<rank4tensor_type> const& fetch(rank4 const matrixType) const
    {
        return rank4tensors.find(matrixType)->second;
    }

    /** Const access to the non-converged scalar variables */
    template <typename... ScalarTps>
    auto fetch(Scalar var0, Scalar var1, ScalarTps... vars) const
    {
        return std::make_tuple(std::cref(scalars.find(var0)->second),
                               std::cref(scalars.find(var1)->second),
                               std::cref(scalars.find(vars)->second)...);
    }

    /** Const access to the non-converged tensor variables */
    template <typename... TensorTps>
    auto fetch(Tensor var0, Tensor var1, TensorTps... vars) const
    {
        return std::make_tuple(std::cref(rank2tensors.find(var0)->second),
                               std::cref(rank2tensors.find(var1)->second),
                               std::cref(rank2tensors.find(vars)->second)...);
    }

    template <typename... rank4_types>
    auto fetch(rank4 const var0, rank4 const var1, rank4_types... vars) const
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

protected:
    // These state variables are committed and reverted depending on the outer
    // simulation loop.  If a nonlinear iteration does not converge then revert
    // the state back to the previous state.  The rank2tensors and scalars fields
    // are the 'unstable' variables and the *Old are the stable variants
    std::unordered_map<Tensor, std::vector<rank2tensor_type>> rank2tensors, rank2tensors_old;
    std::unordered_map<Scalar, std::vector<scalar_type>> scalars, scalars_old;

    std::unordered_map<rank4, std::vector<rank4tensor_type>> rank4tensors;

    std::size_t size;
};
}
