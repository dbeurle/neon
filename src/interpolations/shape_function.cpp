
#include "interpolations/shape_function.hpp"

namespace neon
{
template class shape_function<numerical_quadrature<double>>;
template class shape_function<surface_quadrature>;
template class shape_function<volume_quadrature>;
}
