
#pragma once

#include "InternalVariablesForwards.hpp"

#include "numeric/Tensor.hpp"

#include <functional>
#include <unordered_map>

namespace neon
{
/**
 * InternalVariables stores the internal variables associated with the element
 * quadrature points.  These variables are duplicated and commited to memory
 * when the data is converged to avoid polluting the variable history in the
 * Newton-Raphson method.
 */
template <int spatial_dimension, int voigt_dimension>
class InternalVariables
{
public:
    using tensor = Eigen::Matrix<double, spatial_dimension, spatial_dimension>;
    using matrix = Eigen::Matrix<double, voigt_dimension, voigt_dimension>;

    using Scalars = std::vector<double>;
    using Vectors = std::vector<std::vector<double>>;
    using Tensors = std::vector<tensor>;
    using Matrices = std::vector<matrix>;

public:
    enum class Matrix { TangentOperator };

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
        //  J2Plasticity internal variables
        BackStress,
        KinematicHardening,
        /* Tensors for diffusion applications */
        Conductivity
    };

    enum class Vector { Chains, Segments, ShearModuli, HeatFlux };

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

    /** Implicitly defined move constructor */
    InternalVariables(InternalVariables&&) = default;

    /** Add a number of tensor type variables to the object store */
    template <typename... Variables>
    void add(Tensor const name, Variables... names)
    {
        tensors[name].resize(size, tensor::Zero());
        tensors_old[name].resize(size, tensor::Zero());
        add(names...);
    }

    void add(Tensor const name)
    {
        tensors[name].resize(size, tensor::Zero());
        tensors_old[name].resize(size, tensor::Zero());
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

    /** Allocate matrix internal variable with zeros */
    void add(InternalVariables::Matrix const name) { matrices[name].resize(size, matrix::Zero()); }

    /** Allocate matrix internal variables with provided matrix */
    void add(InternalVariables::Matrix const name, neon::Matrix const initial_matrix)
    {
        matrices[name].resize(size, initial_matrix);
    }

    bool has(Scalar name) const { return scalars.find(name) != scalars.end(); }
    bool has(Tensor name) const { return tensors.find(name) != tensors.end(); }
    bool has(Matrix name) const { return matrices.find(name) != matrices.end(); }

    /** Const access to the converged tensor variables */
    Tensors const& operator[](Tensor const tensorType) const
    {
        return tensors_old.find(tensorType)->second;
    }

    /** Const access to the converged scalar variables */
    Scalars const& operator[](Scalar const scalarType) const
    {
        return scalars_old.find(scalarType)->second;
    }

    /*-------------------------------------------------------------*
     *  Mutable access methods for unconverged internal variables  *
     *-------------------------------------------------------------*/

    /** Mutable access to the non-converged scalar variables */
    Scalars& operator()(Scalar scalarType) { return scalars.find(scalarType)->second; }

    /** Mutable access to the non-converged tensor variables */
    Tensors& operator()(Tensor tensorType) { return tensors.find(tensorType)->second; }

    /** Mutable access to the non-converged matrix variables */
    Matrices& operator()(Matrix matrixType) { return matrices.find(matrixType)->second; }

    /** Mutable access to the non-converged scalar variables */
    template <typename... ScalarTps>
    auto operator()(Scalar var0, Scalar var1, ScalarTps... vars)
    {
        return std::make_tuple(std::ref(scalars.find(var0)->second),
                               std::ref(scalars.find(var1)->second),
                               std::ref(scalars.find(vars)->second)...);
    }

    /** Mutable access to the non-converged tensor variables */
    template <typename... TensorTps>
    auto operator()(Tensor var0, Tensor var1, TensorTps... vars)
    {
        return std::make_tuple(std::ref(tensors.find(var0)->second),
                               std::ref(tensors.find(var1)->second),
                               std::ref(tensors.find(vars)->second)...);
    }

    template <typename... MatrixTps>
    auto operator()(Matrix var0, Matrix var1, MatrixTps... vars)
    {
        return std::make_tuple(std::ref(matrices.find(var0)->second),
                               std::ref(matrices.find(var1)->second),
                               std::ref(matrices.find(vars)->second)...);
    }

    /*-------------------------------------------------------------*
     * Constant access methods for unconverged internal variables  *
     *-------------------------------------------------------------*/

    /** Constant access to the non-converged scalar variables */
    Scalars const& operator()(Scalar scalarType) const { return scalars.find(scalarType)->second; }

    /** Non-mutable access to the non-converged tensor variables */
    Tensors const& operator()(Tensor tensorType) const { return tensors.find(tensorType)->second; }

    /** Non-mutable access to the non-converged matrix variables */
    Matrices const& operator()(Matrix matrixType) const
    {
        return matrices.find(matrixType)->second;
    }

    /** Const access to the non-converged scalar variables */
    template <typename... ScalarTps>
    auto operator()(Scalar var0, Scalar var1, ScalarTps... vars) const
    {
        return std::make_tuple(std::cref(scalars.find(var0)->second),
                               std::cref(scalars.find(var1)->second),
                               std::cref(scalars.find(vars)->second)...);
    }

    /** Const access to the non-converged tensor variables */
    template <typename... TensorTps>
    auto operator()(Tensor var0, Tensor var1, TensorTps... vars) const
    {
        return std::make_tuple(std::cref(tensors.find(var0)->second),
                               std::cref(tensors.find(var1)->second),
                               std::cref(tensors.find(vars)->second)...);
    }

    template <typename... MatrixTps>
    auto operator()(Matrix var0, Matrix var1, MatrixTps... vars) const
    {
        return std::make_tuple(std::cref(matrices.find(var0)->second),
                               std::cref(matrices.find(var1)->second),
                               std::cref(matrices.find(vars)->second)...);
    }

    /** Commit to history when iteration converges */
    void commit()
    {
        tensors_old = tensors;
        scalars_old = scalars;
    }

    /** Revert to the old state when iteration doesn't converge */
    void revert()
    {
        tensors = tensors_old;
        scalars = scalars_old;
    }

protected:
    // These state variables are committed and reverted depending on the outer
    // simulation loop.  If a nonlinear iteration does not converge then revert
    // the state back to the previous state.  The tensors and scalars fields
    // are the 'unstable' variables and the *Old are the stable variants
    std::unordered_map<Tensor, Tensors> tensors, tensors_old;
    std::unordered_map<Scalar, Scalars> scalars, scalars_old;

    std::unordered_map<Matrix, Matrices> matrices;

    std::size_t size;
};
}
