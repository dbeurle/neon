// /*
//  * neon - A finite element solver.
//  *
//  * For licensing please refer to the LICENSE.md file
//  *
//  */
//
// #include "Prism6.hpp"
//
// namespace neon
// {
// Prism6::Prism6(PrismScheme::Rule rule) :
//     integrator(rule, interpolationOrder, dN_dXi, dN_dEta, dN_dZeta)
// {
// 	this->fillShapeFunctionsAndDerivatives();
// }
//
// void Prism6::fillShapeFunctionsAndDerivatives()
// {
//     auto const L = integrator.internalPoints();
//
//     N.resize(6, L);
// 	dN_dXi = dN_dEta = dN_dZeta;
//
// 	for (int l = 0; l < L; ++l)
// 	{
//         auto const tuple = integrator.coordinate(l);
//         auto const r = std::get<0>(tuple);
//         auto const s = std::get<1>(tuple);
//         auto const t = 1.0 - r - s;
//         auto const xi = std::get<2>(tuple);
//
//         N(0, l) = r * (1.0 - xi) / 2.0;
// 		N(1, l) = s * (1.0 - xi) / 2.0;
// 		N(2, l) = t * (1.0 - xi) / 2.0;
// 		N(3, l) = r * (1.0 + xi) / 2.0;
// 		N(4, l) = s * (1.0 + xi) / 2.0;
// 		N(5, l) = t * (1.0 + xi) / 2.0;
//
//         // Compute the shape function derivatives for each quadrature point
// 		dN_dXi(0, l) = (1.0 - xi) / 2.0;
// 		dN_dXi(1, l) = 0.0;
// 		dN_dXi(2, l) = -(1.0 - xi) / 2.0;
// 		dN_dXi(3, l) = (1.0 + xi) / 2.0;
// 		dN_dXi(4, l) = 0.0;
// 		dN_dXi(5, l) = -(1.0 + xi) / 2.0;
//
// 		dN_dEta(0, l) = 0.0;
// 		dN_dEta(1, l) = (1.0 - xi) / 2.0;
// 		dN_dEta(2, l) = -(1.0 - xi) / 2.0;
// 		dN_dEta(3, l) = 0.0;
// 		dN_dEta(4, l) = (1.0 + xi) / 2.0;
// 		dN_dEta(5, l) = -(1.0 + xi) / 2.0;
//
//         // Zeta is the natural coordinate length way on the prism
// 		dN_dZeta(0, l) = -r / 2.0;
// 		dN_dZeta(1, l) = -s / 2.0;
// 		dN_dZeta(2, l) = -t / 2.0;
// 		dN_dZeta(3, l) = r / 2.0;
// 		dN_dZeta(4, l) = s / 2.0;
// 		dN_dZeta(5, l) = t / 2.0;
// 	}
// }
// }
