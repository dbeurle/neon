
#include "modules/SolidMechanicsModule.hpp"

namespace neon
{
SolidMechanicsModule::SolidMechanicsModule(basic_mesh const& mesh,
                                           json const& material,
                                           json const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0]), fem_matrix(fem_mesh, simulation)
{
}
}
