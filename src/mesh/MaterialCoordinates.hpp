
#pragma once

#include "mesh/NodalCoordinates.hpp"

#include "vtkSmartPointer.h"
class vtkPoints;
class vtkDoubleArray;

namespace neon
{
class MaterialCoordinates : public NodalCoordinates
{
public:
    /** Construct this class using a set of initial coordinates */
    MaterialCoordinates(matrix3x const& initial_coordinates);

    /** @return element reference configuration based on the local node numbers*/
    auto const initial_configuration(std::vector<int64> const& nodes) const
    {
        return X(Eigen::placeholders::all, nodes);
    }

    /** @return element current configuration based on the local node numbers*/
    auto const current_configuration(std::vector<int64> const& nodes) const
    {
        return x(Eigen::placeholders::all, nodes);
    }

    /** @param u - displacement vector from initial configuration (x,y,z...) */
    void update_current_configuration(vector const& u);

    /** @param u - displacement vector from initial configuration (x,y...) */
    void update_current_xy_configuration(vector const& u);

    [[nodiscard]] matrix3x displacement() const;

    /** @return a vtk object of the initial coordinates */
    [[nodiscard]] vtkSmartPointer<vtkPoints> vtk_coordinates() const;

    /** @return a vtk array of nodal displacements */
    [[nodiscard]] vtkSmartPointer<vtkDoubleArray> vtk_displacement() const;

protected:
    matrix3x x; //!< Current configuration
};
}
