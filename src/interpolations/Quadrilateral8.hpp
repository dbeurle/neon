// /*
//  * neon - A finite element solver.
//  *
//  * For licensing please refer to the LICENSE.md file
//  *
//  * Copyright Darcy Beurle, 2016.
//  */
//
// #pragma once
//
// #include "SurfaceInterpolation.hpp"
// #include "NumericalIntegration/QuadrilateralScheme.hpp"
//
// namespace neon
// {
// /** A quadrilateral finite element with 8 nodes with an isoparametric formulation */
// class Quadrilateral8: public SurfaceInterpolation
// {
// public:
//
// 	Quadrilateral8(QuadrilateralScheme::Rule rule);
//
// protected:
//
//     /**
//      * Initialize the shape functions of the quadrilateral 8
//      * node element to be the Lagrange polynomials
//      * \f{align*}{
//      * N_1(\xi, \eta) &= \frac{1}{4}\xi\eta(1-\xi)(1-\eta) \\
//      * N_2(\xi, \eta) &= \frac{1}{4}\xi\eta(1+\xi)(1-\eta) \\
//      * N_3(\xi, \eta) &= \frac{1}{4}\xi\eta(1+\xi)(1+\eta) \\
//      * N_4(\xi, \eta) &= \frac{1}{4}\xi\eta(1-\xi)(1+\eta) \\
//      * N_5(\xi, \eta) &= \frac{1}{2}\eta (\xi^2 - 1)(\eta - 1) \\
//      * N_6(\xi, \eta) &= \frac{1}{2}\xi  (\xi + 1  )(1-\eta^2) \\
//      * N_7(\xi, \eta) &= \frac{1}{2}\eta (1 - \xi^2)(\eta + 1) \\
//      * N_8(\xi, \eta) &= \frac{1}{2}\xi  (\xi - 1  )(1 - \eta^2)
//      * \f}
//      */
//     void fillShapeFunctionsAndDerivatives();
//
//     QuadrilateralScheme integrator;
//
// };
// }
