/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 */

#include "Tetrahedron4.hpp"

namespace neon
{
Tetrahedron4::Tetrahedron4() : VolumeInterpolation(nullptr)
{
    dN = Matrix::Zero(4, 3);

    this->setDerivativeMatrix();
    this->initializeSource();
    this->initializeMass();
}

void Tetrahedron4::setDerivativeMatrix()
{
    dN << 1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0,   //
        0.0, 0.0, 1.0,   //
        -1.0, -1.0, -1.0;
}

void Tetrahedron4::initializeSource()
{
    force.setOnes(4);
    force /= 24.0;
}

void Tetrahedron4::initializeMass()
{
    mass.resize(4, 4);
    mass << 2.0, 1.0, 1.0, 1.0, //
        1.0, 2.0, 1.0, 1.0,     //
        1.0, 1.0, 2.0, 1.0,     //
        1.0, 1.0, 1.0, 2.0;
    mass /= 120.0;
}

double Tetrahedron4::computeSixV(Matrix const& nodalCoordinates) const
{
    // Use the determinant to compute the volume of the element
    Eigen::Matrix4d sixV = Eigen::Matrix4d::Ones(4, 4);

    sixV.col(1) = nodalCoordinates.row(0);
    sixV.col(2) = nodalCoordinates.row(1);
    sixV.col(3) = nodalCoordinates.row(2);

    return sixV.determinant();
}

// Matrix Tetrahedron4::computeThermalStiffness( const Matrix& nodalCoordinates,
//                                               const cm::ThermalCEqn& ceqn) const
// {
//     Eigen::Matrix3d Jacobian = nodalCoordinates * dN;
//     Matrix B = (dN * Jacobian.inverse()).transpose();
//     const auto volume = this->computeSixV(nodalCoordinates) / 6.0;
//
//     if (volume < 0)
// 	{
//         Info<< "!Warning: Element with negative determinant detected.\n\tResults may be
//         erroneous";
// 	}
//
//     const auto C = ceqn.getConductivityTensor();
//     return volume * B.transpose() * C * B;
// }
//
// Matrix Tetrahedron4::computeThermalMass( const Matrix& nodalCoordinates,
//                                          const cm::ThermalCEqn& ceqn) const
// {
//     auto const D = ceqn.getSpecificHeat() * ceqn.getDensity();
//     return mass * this->computeSixV(nodalCoordinates) * D;
// }
}
