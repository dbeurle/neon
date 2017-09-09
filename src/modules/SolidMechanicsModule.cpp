
#include "SolidMechanicsModule.hpp"

namespace neon
{
SolidMechanicsModule::SolidMechanicsModule(BasicMesh const& mesh,
                                           Json::Value const& material,
                                           Json::Value const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0]),
      fem_matrix(fem_mesh,
                 Visualisation(simulation["Name"].asString(), fem_mesh, simulation["Visualisation"]),
                 simulation["LinearSolver"],
                 simulation["NonlinearOptions"],
                 simulation["Time"])
{
}
}
