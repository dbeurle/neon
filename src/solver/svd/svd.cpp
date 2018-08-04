
#include "solver/svd/svd.hpp"

namespace neon
{
bdc_svd::bdc_svd(col_matrix const& A) { compute(A); }

void bdc_svd::compute(col_matrix const& A)
{
    A_decomposiion.compute(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    left_vectors = A_decomposiion.matrixU();
    right_vectors = A_decomposiion.matrixV();
    singular_values = A_decomposiion.singularValues();
}

void bdc_svd::compute(col_matrix const& A, int const n)
{
    compute(A);
    left_vectors = left_vectors.leftCols(n).eval();
    right_vectors = A_decomposiion.matrixV().leftCols(n).eval();
    singular_values = A_decomposiion.singularValues().head(n).eval();
}

void bdc_svd::compute(col_matrix const& A, double const tolerance)
{
    compute(A);
    vector normalised_singular_values = singular_values / singular_values(1);
    std::vector<double> singular_values_vector(singular_values.size());
    vector::Map(&singular_values_vector[0],
                normalised_singular_values.size()) = normalised_singular_values;

    auto lower_bound_iterator = std::lower_bound(singular_values_vector.begin(),
                                                 singular_values_vector.end(),
                                                 tolerance,
                                                 [](double const& x, double const d) {
                                                     return x >= d;
                                                 });
    int n = std::distance(singular_values_vector.begin(), lower_bound_iterator);
    left_vectors = A_decomposiion.matrixU().leftCols(n).eval();
    right_vectors = A_decomposiion.matrixV().leftCols(n).eval();
    singular_values = A_decomposiion.singularValues().head(n).eval();
}

col_matrix const& bdc_svd::left() const noexcept { return left_vectors; }

col_matrix const& bdc_svd::right() const noexcept { return right_vectors; }

vector const& bdc_svd::values() const noexcept { return singular_values; }

void bdc_svd::solve(vector& x, vector& b) const noexcept { x = A_decomposiion.solve(b); }
}
