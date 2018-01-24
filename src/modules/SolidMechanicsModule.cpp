
#include "SolidMechanicsModule.hpp"

namespace neon
{
SolidMechanicsModule::SolidMechanicsModule(BasicMesh const& mesh,
                                           json const& material,
                                           json const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0]), fem_matrix(fem_mesh, simulation)
{
}

PlaneMechanicsModule::PlaneMechanicsModule(BasicMesh const& mesh,
                                           json const& material,
                                           json const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0]), fem_matrix(fem_mesh, simulation)
{
}
}
