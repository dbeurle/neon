
#include "solver/svd/svd.hpp"

namespace neon
{
bdc_svd::bdc_svd(col_matrix const& A) { compute(A); }

void bdc_svd::compute(col_matrix const& A)
{
    decomposition.compute(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    left_vectors = decomposition.matrixU();
    right_vectors = decomposition.matrixV();
    singular_values = decomposition.singularValues();
}

void bdc_svd::compute(col_matrix const& A, std::int64_t const n)
{
    compute(A);
    left_vectors = left_vectors.leftCols(n).eval();
    right_vectors = decomposition.matrixV().leftCols(n).eval();
    singular_values = decomposition.singularValues().head(n).eval();
}

void bdc_svd::compute(col_matrix const& A, double const tolerance)
{
    compute(A);

    vector normalised_singular_values = singular_values / singular_values(1);
    std::vector<double> singular_values_vector(singular_values.size());

    vector::Map(singular_values_vector.data(),
                normalised_singular_values.size()) = normalised_singular_values;

    auto const n = std::distance(begin(singular_values_vector),
                                 std::lower_bound(begin(singular_values_vector),
                                                  singular_values_vector.end(),
                                                  tolerance,
                                                  std::greater_equal{}));

    left_vectors = decomposition.matrixU().leftCols(n).eval();
    right_vectors = decomposition.matrixV().leftCols(n).eval();
    singular_values = decomposition.singularValues().head(n).eval();
}

col_matrix const& bdc_svd::left() const noexcept { return left_vectors; }

col_matrix const& bdc_svd::right() const noexcept { return right_vectors; }

vector const& bdc_svd::values() const noexcept { return singular_values; }

void bdc_svd::solve(vector& x, vector const& b) const noexcept { x = decomposition.solve(b); }
}
