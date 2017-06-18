
#include "NumericalQuadrature.hpp"

namespace neon
{
/** UnitSphereQuadrature implements the symmetric quadrature rules */
class UnitSphereQuadrature : public VolumeQuadrature
{
public:
    /** Fill the quadrature coordinates and weightings */
    UnitSphereQuadrature();
};
}
