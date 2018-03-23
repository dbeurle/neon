
#pragma once

#include "shape_function.hpp"
#include "quadrature/tetrahedron_quadrature.hpp"

namespace neon
{
/// tetrahedron4 is an isoparametric simplex 4 node element numerical integration.
class tetrahedron4 : public volume_interpolation
{
public:
    tetrahedron4(tetrahedron_quadrature::Rule rule = tetrahedron_quadrature::Rule::OnePoint);

    virtual int nodes() const override final { return 4; }

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t) &= r \\
     * N_2(r, s, t) &= s \\
     * N_3(r, s, t) &= t \\
     * N_4(r, s, t) &= 1 - r - s - t
     * \f}
     */
    void precompute_shape_functions();
};

/** Isoparametric quadratric tetrahedral element (10 nodes) */
class tetrahedron10 : public volume_interpolation
{
public:
    tetrahedron10(tetrahedron_quadrature::Rule rule);

    virtual int nodes() const override final { return 10; }

protected:
    /**
     * Initialize the shape functions to be given the following trial functions
     * using tetrahedral coordinates:
     * \f{align*}{
     * N_1(r,s,t,u) &= r(2r-1) \\
     * N_2(r,s,t,u) &= s(2s-1) \\
     * N_3(r,s,t,u) &= t(2t-1) \\
     * N_4(r,s,t,u) &= u(2u-1) \\
     * N_5(r,s,t,u) &= 4rs \\
     * N_6(r,s,t,u) &= 4st \\
     * N_7(r,s,t,u) &= 4tu \\
     * N_8(r,s,t,u) &= 4ru \\
     * N_9(r,s,t,u) &= 4rt \\
     * N_{10}(r,s,t,u) &= 4su
     * \f}
     * And the derivative of the shape functions in the natural coordinates:
     * \f{align*}{
     * \frac{\partial{N_1}}{\partial{r}} &= 4r-1        &
     * \frac{\partial{N_1}}{\partial{s}}    &= 0
     * & \frac{\partial{N_1}}{\partial{t}}    &= 0           \\
     * \frac{\partial{N_2}}{\partial{r}} &= 0           &
     * \frac{\partial{N_2}}{\partial{s}}    &=
     * 4s-1        & \frac{\partial{N_2}}{\partial{t}}    &= 0           \\
     * \frac{\partial{N_3}}{\partial{r}} &= 0           &
     * \frac{\partial{N_3}}{\partial{s}}    &= 0
     * & \frac{\partial{N_3}}{\partial{t}}    &= 4t-1        \\
     * \frac{\partial{N_4}}{\partial{r}} &= -3+4r+4s+4t &
     * \frac{\partial{N_4}}{\partial{s}}    &=
     * -3+4r+4s+4t & \frac{\partial{N_4}}{\partial{t}}    &= -3+4r+4s+4t \\
     * \frac{\partial{N_5}}{\partial{r}} &= 0           &
     * \frac{\partial{N_5}}{\partial{s}}    &= 4r
     * & \frac{\partial{N_5}}{\partial{t}}    &= 0           \\
     * \frac{\partial{N_6}}{\partial{r}} &= 4s          &
     * \frac{\partial{N_6}}{\partial{s}}    &= 4t
     * & \frac{\partial{N_6}}{\partial{t}}    &= 4s          \\
     * \frac{\partial{N_7}}{\partial{r}} &= -4t         &
     * \frac{\partial{N_7}}{\partial{s}}    &=
     * -4t         & \frac{\partial{N_7}}{\partial{t}}    &= 4-4r-4s-8t  \\
     * \frac{\partial{N_8}}{\partial{r}} &= 4-4s-8r-4t  &
     * \frac{\partial{N_8}}{\partial{s}}    &=
     * -4r         & \frac{\partial{N_8}}{\partial{t}}    &= -4r         \\
     * \frac{\partial{N_9}}{\partial{r}} &= 4t          &
     * \frac{\partial{N_9}}{\partial{s}}    &= 0
     * & \frac{\partial{N_9}}{\partial{t}}    &= 4r          \\
     * \frac{\partial{N_{10}}}{\partial{r}} &= -4s      &
     * \frac{\partial{N_{10}}}{\partial{s}} &=
     * 4-4r-8s-4t  & \frac{\partial{N_{10}}}{\partial{t}} &= -4s
     * \f}
     */
    void precompute_shape_functions();
};
}
