
#include "solver/svd/svd.hpp"

namespace neon
{
randomised_svd::randomised_svd(col_matrix const& A) { compute(A); }

void randomised_svd::compute(col_matrix const& A)
{
    A_decomposiion.compute(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    left_vectors = A_decomposiion.matrixU();
    right_vectors = A_decomposiion.matrixV();
    singular_values = A_decomposiion.singularValues();
}

void randomised_svd::compute(col_matrix const& A, int const n)
{
    int cols = A.cols();
    int random_matrix_cols = std::min(2 * n, cols);

    col_matrix random_matrix = col_matrix::Random(cols, random_matrix_cols);
    col_matrix range_approximation = A * random_matrix;

    compute(range_approximation); // TODO could we comput only left_vectors?
    col_matrix A_projection = left_vectors.transpose() * A;
    col_matrix orthonormal_basis = left_vectors;

    compute(A_projection);
    left_vectors = orthonormal_basis * left_vectors.eval();

    int k = std::min(n, int(left_vectors.cols()));
    left_vectors = left_vectors.leftCols(k).eval();
    right_vectors = A_decomposiion.matrixV().leftCols(k).eval();
    singular_values = A_decomposiion.singularValues().head(k).eval();
}

void randomised_svd::compute(col_matrix const& A, double const tolerance)
{
    compute(A, int(A.cols()));
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

col_matrix const& randomised_svd::left() const noexcept { return left_vectors; }

col_matrix const& randomised_svd::right() const noexcept { return right_vectors; }

vector const& randomised_svd::values() const noexcept { return singular_values; }

void randomised_svd::solve(vector& x, vector& b) const noexcept { x = A_decomposiion.solve(b); }
}
