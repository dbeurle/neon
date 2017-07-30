
#pragma once

#include <functional>
#include <unordered_map>
#include <vector>

#include "numeric/DenseTypes.hpp"
#include "numeric/Tensor.hpp"

namespace neon
{
/**
 * InternalVariables stores the internal variables associated with the element
 * quadrature points.  These variables are duplicated and commited to memory
 * when the data is converged to avoid polluting the variable history in the
 * Newton-Raphson method.
 */
class InternalVariables
{
public:
    using Scalars = std::vector<double>;
    using Tensors = std::vector<Matrix3>;
    using Matrices = std::vector<CMatrix>;

    enum class Matrix { TruesdellModuli };

    enum class Tensor {
        // Stress measures
        Cauchy,
        Kirchhoff,
        PiolaKirchhoff1,
        PiolaKirchhoff2,
        // Deformation measures
        LinearisedStrain,
        LinearisedPlasticStrain, // Small strain
        RateOfDeformation,
        RateOfDeformationPlastic,
        DeformationGradient,
        DeformationGradientPlastic,
        DisplacementGradient,
        GreenLagrange
    };

    enum class Scalar {
        Chains,   // Chains for the micromechanical model
        Segments, // Segments for the micromechanical model
        VonMisesStress,
        EffectivePlasticStrain,
        DetF0, // Reference Jacobian determinant
        DetF   // Updated Jacobian determinant
    };

public:
    InternalVariables(std::size_t const size) : size(size) {}

    /** Delete copy constructor to prevent references moving */
    InternalVariables(InternalVariables const&) = delete;

    /** Implicitly defined move constructor */
    InternalVariables(InternalVariables&&) = default;

    /** Add a number of tensor type variables to the object store */
    template <typename... Variables>
    void add(Tensor name, Variables... names);

    void add(Tensor name);

    /** Add a number of scalar type variables to the object store */
    template <typename... Variables>
    void add(Scalar name, Variables... names);

    void add(Scalar name);

    /**
     * @param size Size of internal buffer
     * @param name Name of variables
     * @param rowcol Number of rows (or columns) in the square matrix
     */
    void add(Matrix name, int rowcol);

    /**
     * @param size Size of internal buffer
     * @param name Name of variables
     * @param init Initial matrix
     */
    void add(Matrix name, neon::Matrix init);

    bool has(Scalar name) const { return scalars.find(name) != scalars.end(); }
    bool has(Tensor name) const { return tensors.find(name) != tensors.end(); }
    bool has(Matrix name) const { return matrices.find(name) != matrices.end(); }

    /** Const access to the converged tensor variables */
    Tensors const& operator[](Tensor tensorType) const;

    /** Const access to the converged scalar variables */
    Scalars const& operator[](Scalar scalarType) const;

    /** Const access to the converged matrix variables */
    Matrices const& operator[](Matrix matrixType) const;

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
    Scalars const& operator()(Scalar scalarType) const
    {
        return scalars.find(scalarType)->second;
    }

    /** Non-mutable access to the non-converged tensor variables */
    Tensors const& operator()(Tensor tensorType) const
    {
        return tensors.find(tensorType)->second;
    }

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
    void commit();

    /** Revert to the old state when iteration doesn't converge */
    void revert();

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

// Allocation methods

template <typename... Variables>
inline void InternalVariables::add(InternalVariables::Tensor name, Variables... names)
{
    tensors[name].resize(size, Matrix3::Zero());
    tensors_old[name].resize(size, Matrix3::Zero());
    add(names...);
}

inline void InternalVariables::add(InternalVariables::Tensor name)
{
    tensors[name].resize(size, Matrix3::Zero());
    tensors_old[name].resize(size, Matrix3::Zero());
}

template <typename... Variables>
inline void InternalVariables::add(InternalVariables::Scalar name, Variables... names)
{
    scalars[name].resize(size, 0.0);
    scalars_old[name].resize(size, 0.0);
    add(names...);
}

inline void InternalVariables::add(InternalVariables::Scalar name)
{
    scalars[name].resize(size, 0.0);
    scalars_old[name].resize(size, 0.0);
}

inline void InternalVariables::add(InternalVariables::Matrix name, int rowcol)
{
    matrices[name].resize(size, neon::Matrix::Zero(rowcol, rowcol));
}

inline void InternalVariables::add(InternalVariables::Matrix name, neon::Matrix init)
{
    matrices[name].resize(size, init);
}

// Converged results

inline InternalVariables::Tensors const& InternalVariables::operator[](
    InternalVariables::Tensor tensorType) const
{
    return tensors_old.find(tensorType)->second;
}

inline InternalVariables::Scalars const& InternalVariables::operator[](Scalar scalarType) const
{
    return scalars_old.find(scalarType)->second;
}

// Version control of internal state variables

inline void InternalVariables::commit()
{
    tensors_old = tensors;
    scalars_old = scalars;
}

inline void InternalVariables::revert()
{
    tensors = tensors_old;
    scalars = scalars_old;
}
}
