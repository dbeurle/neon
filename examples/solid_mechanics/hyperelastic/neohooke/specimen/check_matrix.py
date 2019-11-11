import scipy.io
import scipy.sparse.csgraph
from scipy.sparse import csr_matrix

import numpy
import matplotlib.pyplot as plt

A = scipy.io.mmread("sparse_matrix.mtx").toarray()
b = scipy.io.mmread("right_hand_side.mtx")

delta_u = numpy.linalg.solve(A, b)

print(numpy.linalg.norm(delta_u))

permutation = scipy.sparse.csgraph.reverse_cuthill_mckee(csr_matrix(A))

reordered_A = numpy.squeeze(numpy.array(A))
reordered_A = reordered_A[permutation, :]
reordered_A = reordered_A[:, permutation]

(w, v) = numpy.linalg.eig(A)

w = numpy.sort(w)

print("condition number", w[-1] / w[0])

plt.figure(1)
plt.spy(A)

plt.figure(2)
plt.spy(reordered_A)

plt.show()
