//
// #pragma once
//
// #include "mesh/common/Neumann.hpp"
//
// #include "interpolations/SurfaceInterpolation.hpp"
// #include "mesh/NodalCoordinates.hpp"
// #include "mesh/diffusion/Submesh.hpp"
//
// namespace neon::diffusion
// {
// /**
//  * SurfaceFlux is a group of the same elements on the same boundary. These
//  * elements are responsible for computing their elemental right hand side contributions
//  * with the corresponding shape function.  These are required to be stored in a parent
//  * container with the other groups from the collective boundary \sa SurfaceFluxBoundary
//  */
// using SurfaceFlux = SurfaceLoad<SurfaceInterpolation>;
//
// using SurfaceFluxBoundary = SurfaceLoadBoundary<1, SurfaceFlux>;
// }
