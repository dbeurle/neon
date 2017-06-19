// /*
//  * neon - A finite element solver.
//  *
//  * For licensing please refer to the LICENSE.md file
//  *
//  */
//
// #include "Quadrilateral8.hpp"
//
// namespace neon
// {
// Quadrilateral8::Quadrilateral8(QuadrilateralScheme::Rule rule):
//     integrator(rule, 2, dN_dXi, dN_dEta)
// {
// 	this->fillShapeFunctionsAndDerivatives();
// }
//
// void Quadrilateral8::fillShapeFunctionsAndDerivatives()
// {
//     auto const abscissas = integrator.internalPoints();
//
//     N.resize(8, abscissas);
//
//     dN_dXi = dN_dEta = N;
//
// 	for (int l = 0; l < abscissas; ++l)
// 	{
//         auto const naturalCoordinates = integrator.coordinate(l);
//
//         auto const xi  = std::get<0>(naturalCoordinates);
//         auto const eta = std::get<1>(naturalCoordinates);
//
//         N(0, l) = 0.25 * xi * eta * (xi - 1.0) * (eta - 1.0);
//         N(1, l) = 0.25 * xi * eta * (xi + 1.0) * (eta - 1.0);
//         N(2, l) = 0.25 * xi * eta * (xi + 1.0) * (eta + 1.0);
//         N(3, l) = 0.25 * xi * eta * (xi - 1.0) * (eta + 1.0);
//         N(4, l) = 0.5 * eta * (1.0 - xi * xi)  * (eta - 1.0);
//         N(5, l) = 0.5 * xi  * (xi + 1.0) * (1.0 - eta * eta);
//         N(6, l) = 0.5 * eta * (1.0 - xi * xi)  * (eta + 1.0);
//         N(7, l) = 0.5 * xi  * (xi - 1.0) * (1.0 - eta * eta);
//
//         dN_dXi(0, l) = 0.25 * eta * (2.0 * xi + 1.0) * (eta - 1.0);
//         dN_dXi(1, l) = 0.25 * eta * (2.0 * xi - 1.0) * (eta - 1.0);
//         dN_dXi(2, l) = 0.25 * eta * (2.0 * xi + 1.0) * (eta + 1.0);
//         dN_dXi(3, l) = 0.25 * eta * (2.0 * xi - 1.0) * (eta + 1.0);
//         dN_dXi(4, l) = -xi * eta * (eta - 1.0);
//         dN_dXi(5, l) = 0.5 * (2.0 * xi + 1.0) * (1.0 - eta * eta);
//         dN_dXi(6, l) = -xi * eta * (eta + 1.0);
//         dN_dXi(7, l) = 0.5 * (2.0 * xi - 1.0) * (1.0 - eta * eta);
//
//         dN_dEta(0, l) = 0.25 * xi * (xi - 1.0) * (2.0 * eta - 1.0);
//         dN_dEta(1, l) = 0.25 * xi * (xi + 1.0) * (2.0 * eta - 1.0);
//         dN_dEta(2, l) = 0.25 * xi * (xi + 1.0) * (2.0 * eta + 1.0);
//         dN_dEta(3, l) = 0.25 * xi * (xi - 1.0) * (2.0 * eta + 1.0);
//         dN_dEta(4, l) = 0.5 * (1.0 - xi * xi) * (2.0 * eta - 1.0);
//         dN_dEta(5, l) = -xi * eta * (xi + 1.0);
//         dN_dEta(6, l) = 0.5 * (1.0 - xi * xi) * (2.0 * eta + 1.0);
//         dN_dEta(7, l) = -xi * eta * (eta - 1.0);
// 	}
// }
// }
