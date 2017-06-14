
#pragma once

#include "DenseTypes.hpp"
#include <iomanip>
#include <iostream>

namespace neon
{
/** TensorExpression is a CRTP base class for Tensors */
template <typename T>
class TensorExpression
{
public:
    operator T&() { return static_cast<T&>(*this); }

    operator T const&() const { return static_cast<const T&>(*this); }

    double operator()(int i) const { return static_cast<T const&>(*this)(i); }
};

/**
 * Tensor is an expression templated symmetric tensor class for use with
 * stress and strain tensors that possess symmetry.  The available methods
 * reflect tensor operations naturally and facilitate the computation of
 * stress update algorithms which follow the mathematical formulation.
 */
class Tensor : public TensorExpression<Tensor>
{
public:
    /** Tnitialize tensor to zero */
    Tensor() : tau(Matrix3::Zero()) {}

    Tensor(Matrix3 const& tau) : tau(tau) {}

    Tensor(decltype(Matrix3::Identity()) I) : tau(I) {}

    template <typename LeftProduct_Tp, typename RightProduct_Tp>
    Tensor(Eigen::Product<LeftProduct_Tp, RightProduct_Tp> const&& rvalue)
    {
        tau(0, 0) = rvalue(0);
        tau(1, 1) = rvalue(1);
        tau(2, 2) = rvalue(2);
    }

    /** Construct the tensor from an Eigen vector in Voigt form */
    Tensor(Vector const& v) { std::cout << "Constructing from vector"; }

    /** Copy constructor for expression templating magic */
    template <typename Tp>
    Tensor(TensorExpression<Tp> const& tensor)
    {
        for (int i = 0; i < 6; ++i) sigma[i] = tensor(i);
    }

    /** Assigment operator */
    template <typename Tp>
    void operator=(Tp rvalue)
    {
        tau = rvalue;
    }

    template <typename Tp>
    void operator+=(TensorExpression<Tp> const& tensor)
    {
        for (int i = 0; i < 6; ++i) sigma[i] += tensor(i);
    }

    template <typename Tp>
    void operator-=(TensorExpression<Tp> const& tensor)
    {
        for (int i = 0; i < 6; ++i) sigma[i] -= tensor(i);
    }

    double operator()(int i) const { return sigma[i]; }

    Matrix3 const& matrix() const { return tau; }

    Vector voigt() const
    {
        Vector v(6);
        v << sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5];
        return v;
    }

    void voigt(const Vector& t)
    {
        sigma[0] = t(0);
        sigma[1] = t(1);
        sigma[2] = t(2);
        sigma[3] = t(3);
        sigma[4] = t(4);
        sigma[5] = t(5);
    }

    void voigt(double v[6]) { std::copy_n(v, 6, sigma); }

    /** Return the coefficient wise square root */
    // Tensor cwiseSqrt() const { return Tensor(tau.cwiseSqrt()); }

    /**
     * Tensor deviatoric() is the deviatoric stress tensor
     * given by subtraction of the hydrostatic stress tensor
     * from the Cauchy stress tensor:
     * \f{align*}{
     * s_{ij} &= \sigma_{ij} - \frac{\sigma_{kk}}{3} \delta_{ij}\\
     * &= \begin{bmatrix}\sigma_{11} - \pi & \sigma_{12} & \sigma_{13} \\
     * 					 \sigma_{21} & \sigma_{22} - \pi & \sigma_{23} \\
     * 					 \sigma_{31} & \sigma_{32} & \sigma_{33} - \pi \\
     * 					 \end{bmatrix}
     * \f}
     * @return deviatoric part of this Tensor
     */
    Tensor deviatoric() const
    {
        Matrix3 ttau = tau;
        ttau -= (Eigen::Vector3d::Ones() * tau.trace() / 3.0).asDiagonal();
        return Tensor(ttau);
    }

    void print() const { std::cout << "\n" << tau << "\n"; }

    /**
     * trace() calculates the trace of this Tensor
     * \f{align*}{tr(A) = \sum\limits_{i=1}^n a_{ii} \f}
     */
    double trace() const { return tau.trace(); }

    /**
     * determinant() calculates the determinant of this Tensor
     * \f{align*}{\det(\boldsymbol{\sigma}) &= \sigma_{11} \sigma_{22} \sigma_{33}
     * + 2 \sigma_{12} \sigma_{23} \sigma_{31} - \sigma_{12}^2 \sigma_{33}
     * - \sigma_{23}^2 \sigma_{11} - \sigma_{31}^2 \sigma_{22} \f}
     * @return determinant of this Tensor
     */
    double determinant() const { return tau.determinant(); }

    /**
     * I1 returns the coefficient I1, the first stress invariant,
     * which is equal to the trace
     * @return First invariant
     */
    double I1() const { return this->trace(); }

    /**
     * I2 returns the coefficient I2, the second stress invariant,
     * which is calculated by:
     * \f{align*}{
     * I_2 &= \frac{1}{2} \left( (tr \tau)^2 - tr(\tau \tau) \right)
     * \f}
     * @return Second invariant
     */
    double I2() const { return 0.5 * (std::pow(tau.trace(), 2) - (tau * tau).trace()); }

    /** @return Third invariant */
    double I3() const { return this->determinant(); }

    /** Compute the norm of the symmetric tensor */
    double norm() const { return tau.norm(); }

    /** ddot performs the colon operator on itself and returns a scalar value */
    double ddot() const { return (tau.array() * tau.array()).sum(); }

protected:
    alignas(16) double sigma[6];
    Matrix3 tau;
};

// Expression templating engines
template <typename T1, typename T2>
class TensorDifference : public TensorExpression<TensorDifference<T1, T2>>
{
public:
    TensorDifference(TensorExpression<T1> const& t1, TensorExpression<T2> const& t2)
        : tensor1(t1), tensor2(t2)
    {
    }

    double operator()(int i) const { return tensor1(i) - tensor2(i); }

private:
    T1 const& tensor1;
    T2 const& tensor2;
};

template <typename T1, typename T2>
class TensorSum : public TensorExpression<TensorSum<T1, T2>>
{
public:
    TensorSum(TensorExpression<T1> const& t1, TensorExpression<T2> const& t2)
        : tensor1(t1), tensor2(t2)
    {
    }

    double operator()(int i) const { return tensor1(i) + tensor2(i); }

private:
    T1 const& tensor1;
    T2 const& tensor2;
};

template <typename T>
class TensorScaled : public TensorExpression<TensorScaled<T>>
{
public:
    TensorScaled(double a, TensorExpression<T> const& t) : scale(a), tensor(t) {}

    double operator()(int i) const { return scale * tensor(i); }

private:
    double scale;
    T const& tensor;
};

template <typename T1, typename T2>
class TensorDDot : public TensorExpression<TensorDDot<T1, T2>>
{
public:
    TensorDDot(TensorExpression<T1> const& s1, TensorExpression<T2> const& s2) : s1(s1), s2(s2) {}

    double operator()(int i) const { return s1(i) * s2(i); }

private:
    T1 const& s1;
    T2 const& s2;
};

template <typename T>
class TensorQuotient : public TensorExpression<TensorQuotient<T>>
{
public:
    TensorQuotient(TensorExpression<T> const& t, double a) : tensor(t), quotient(a) {}

    double operator()(int i) const { return tensor(i) / quotient; }

private:
    T const& tensor;
    double quotient;
};

// Overloaded operators and helper functions

template <typename T1, typename T2>
TensorDDot<T1, T2> const ddot(TensorExpression<T1> const& t1, TensorExpression<T2> const& t2)
{
    return TensorDDot<T1, T2>(t1, t2);
}

template <typename T1, typename T2>
TensorDifference<T1, T2> const operator-(TensorExpression<T1> const& t1,
                                         TensorExpression<T2> const& t2)
{
    return TensorDifference<T1, T2>(t1, t2);
}

template <typename T1, typename T2>
TensorSum<T1, T2> const operator+(TensorExpression<T1> const& t1, TensorExpression<T2> const& t2)
{
    return TensorSum<T1, T2>(t1, t2);
}

template <typename T>
TensorScaled<T> const operator*(double alpha, TensorExpression<T> const& t)
{
    return TensorScaled<T>(alpha, t);
}

template <typename T>
TensorScaled<T> const operator*(TensorExpression<T> const& t, double alpha)
{
    return TensorScaled<T>(alpha, t);
}

template <typename T>
TensorQuotient<T> const operator/(TensorExpression<T> const& t, double alpha)
{
    return TensorQuotient<T>(t, alpha);
}
}
