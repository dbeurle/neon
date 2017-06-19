// /*
//  * neon - A finite element solver.
//  *
//  * For licensing please refer to the LICENSE.md file
//  *
//  */
//
// #pragma once
//
// #include "VolumeInterpolation.hpp"
// #include "NumericalIntegration/PrismQuadrature.hpp"
//
// namespace neon
// {
// /** Prism6 finite element with 6 nodal points and isoparametric formulation */
// class Prism6: public VolumeInterpolation
// {
// public:
//
//     static constexpr auto interpolationOrder = 1;
//
// 	Prism6(PrismScheme::Rule rule);
//
// protected:
//
//     /**
//      * Initialize the shape functions to be given the following interpolation
//      * functions using triangular and Cartesian coordinates:
//      * \f{align*}{
//      * N_1(r,s,t,\xi) &= \frac{r(1 - \xi)}{2} \\
//      * N_2(r,s,t,\xi) &= \frac{s(1 - \xi)}{2} \\
//      * N_3(r,s,t,\xi) &= \frac{t(1 - \xi)}{2} \\
//      * N_4(r,s,t,\xi) &= \frac{r(1 + \xi)}{2} \\
//      * N_5(r,s,t,\xi) &= \frac{s(1 + \xi)}{2} \\
//      * N_6(r,s,t,\xi) &= \frac{t(1 + \xi)}{2}
//      * \f}
//      * The derivatives are given by:
//      * \f{align*}{
//      * \frac{\partial{N_1}}{\partial{r}} &= \frac{1}{2}(1-\xi)  & \frac{\partial{N_1}}{\partial{s}}    &= 0                   & \frac{\partial{N_1}}{\partial{\xi}}    &= -\frac{1}{2}r \\
//      * \frac{\partial{N_2}}{\partial{r}} &= 0                   & \frac{\partial{N_2}}{\partial{s}}    &= \frac{1}{2}(1-\xi)  & \frac{\partial{N_2}}{\partial{\xi}}    &= -\frac{1}{2}s \\
//      * \frac{\partial{N_3}}{\partial{r}} &= -\frac{1}{2}(1-\xi) & \frac{\partial{N_3}}{\partial{s}}    &= -\frac{1}{2}(1-\xi) & \frac{\partial{N_3}}{\partial{\xi}}    &= -\frac{1}{2}t \\
//      * \frac{\partial{N_4}}{\partial{r}} &= \frac{1}{2}(1+\xi)  & \frac{\partial{N_4}}{\partial{s}}    &= 0                   & \frac{\partial{N_4}}{\partial{\xi}}    &= \frac{1}{2}r  \\
//      * \frac{\partial{N_5}}{\partial{r}} &= 0                   & \frac{\partial{N_5}}{\partial{s}}    &= \frac{1}{2}(1+\xi)  & \frac{\partial{N_5}}{\partial{\xi}}    &= \frac{1}{2}s  \\
//      * \frac{\partial{N_6}}{\partial{r}} &= -\frac{1}{2}(1+\xi) & \frac{\partial{N_6}}{\partial{s}}    &= -\frac{1}{2}(1+\xi) & \frac{\partial{N_6}}{\partial{\xi}}    &= \frac{1}{2}t  \\
//      * \f}
//      * TODO Write the rest of these out
//      */
// 	void fillShapeFunctionsAndDerivatives();
//
//     PrismScheme integrator;
//
// };
// }
