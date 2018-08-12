
#include "solver/svd/svd.hpp"

namespace neon
{
randomised_svd::randomised_svd(col_matrix const& A) { compute(A); }

void randomised_svd::compute(col_matrix const& A)
{
    decomposition.compute(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    left_vectors = decomposition.matrixU();
    right_vectors = decomposition.matrixV();
    singular_values = decomposition.singularValues();
}

void randomised_svd::compute(col_matrix const& A, std::int64_t const n)
{
    auto const cols = A.cols();
    auto const random_matrix_cols = std::min(2 * n, cols);

    col_matrix random_matrix = col_matrix::Random(cols, random_matrix_cols);
    col_matrix range_approximation = A * random_matrix;

    compute(range_approximation); // TODO could we compute only left_vectors?
    col_matrix A_projection = left_vectors.transpose() * A;
    col_matrix orthonormal_basis = left_vectors;

    compute(A_projection);
    left_vectors = orthonormal_basis * left_vectors.eval();

    auto const k = std::min(n, left_vectors.cols());

    left_vectors = left_vectors.leftCols(k).eval();
    right_vectors = decomposition.matrixV().leftCols(k).eval();
    singular_values = decomposition.singularValues().head(k).eval();
}

void randomised_svd::compute(col_matrix const& A, double const tolerance)
{
    compute(A, A.cols());

    vector normalised_singular_values = singular_values / singular_values(1);
    std::vector<double> singular_values_vector(singular_values.size());

    vector::Map(singular_values_vector.data(),
                normalised_singular_values.size()) = normalised_singular_values;

    auto const n = std::distance(begin(singular_values_vector),
                                 std::lower_bound(begin(singular_values_vector),
                                                  end(singular_values_vector),
                                                  tolerance,
                                                  std::greater_equal{}));

    left_vectors = decomposition.matrixU().leftCols(n).eval();
    right_vectors = decomposition.matrixV().leftCols(n).eval();
    singular_values = decomposition.singularValues().head(n).eval();
}

col_matrix const& randomised_svd::left() const noexcept { return left_vectors; }

col_matrix const& randomised_svd::right() const noexcept { return right_vectors; }

vector const& randomised_svd::values() const noexcept { return singular_values; }

void randomised_svd::solve(vector& x, vector const& b) const noexcept
{
    x = decomposition.solve(b);
}
}
