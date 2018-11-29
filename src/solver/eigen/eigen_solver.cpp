
#include "eigen_solver.hpp"

namespace neon
{
eigen_solver::eigen_solver(std::int64_t const values_to_extract, eigen_spectrum const spectrum)
    : values_to_extract{values_to_extract}, m_spectrum{spectrum}
{
}
}
