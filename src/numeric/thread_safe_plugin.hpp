
/** Add a scalar to the value of the matrix at position \a i, \a j in a thread safe fashion.
 *
 * It is assumed that the element exists, otherwise an exception std::domain_error is thrown.
 *
 * This is a O(log(nnz_j)) operation (binary search) with the lock placed around the update
 * of the coefficient matrix
 */
void coefficient_update(Index const row, Index const col, Scalar const value)
{
    eigen_assert(row >= 0 && row < rows() && col >= 0 && col < cols());

    Index const outer = IsRowMajor ? row : col;
    Index const inner = IsRowMajor ? col : row;

    Index const start = m_outerIndex[outer];
    Index const end = m_outerIndex[outer + 1];

    eigen_assert(end >= start && "you probably called coefficient_update on a non finalized matrix");

    Index const p = m_data.searchLowerIndex(start, end - 1, StorageIndex(inner));

// std::lock_guard<std::mutex> lock(coefficent_update_mutex);
#pragma omp atomic
    m_data.value(p) += value;
    // }
}
