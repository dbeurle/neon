
from sympy import *
import numpy as np

from sympy.tensor.array import Array

init_printing()

# Compute the tangent stiffness matrix for J2 plasticity using
# the Truesdell rate of the Cauchy stress.

# This script symbolically computes the factors for the 2 and 3D cases
# and is used in the hypoelastic-plastic models in the neon finite element code.

# Compute:  Cdash = Cstar + sigma *outer_product* I
# Result: Tangent stiffness in Voigt notation


#--------------------------------------------------------------#
#               Two dimensional case (plane strain)            #
#--------------------------------------------------------------#

# Store the mapping for the Voigt notation from Belytschko
voigt_map = [[0, 0], [1, 1], [0, 1]]

# C_dash
σ = MatrixSymbol('σ', 2, 2)

pprint(σ.as_explicit())

Cdash = zeros(3, 3)

for a in range(0, 3):
    i, j = voigt_map[a]
    for b in range(0, 3):
        k, l = voigt_map[b]

        δik = int(i == k)
        δjl = int(j == l)
        δil = int(i == l)
        δjk = int(j == k)

        Cdash[a, b] = Rational(1, 2) * (δik*σ[j,l] + δil*σ[j,k] + δjk*σ[i,l] + δjl*σ[i,k])

print('Plane strain case\n')
pprint(Cdash)

#--------------------------------------------------------------#
#                    Three dimensional case                    #
#--------------------------------------------------------------#

# Store the mapping for the Voigt notation from Belytschko
voigt_map = [[0, 0], [1, 1], [2, 2], [1, 2], [0, 2], [0, 1]]

# C_dash
σ = MatrixSymbol('σ', 3, 3)

Cdash = zeros(6, 6)

for a in range(0, 6):
    i, j = voigt_map[a]
    for b in range(0, 6):
        k, l = voigt_map[b]

        δik = S(int(i == k))
        δjl = S(int(j == l))
        δil = S(int(i == l))
        δjk = S(int(j == k))

        Cdash[a, b] = Rational(1, 2) * (δik*σ[j,l] + δil*σ[j,k] + δjk*σ[i,l] + δjl*σ[i,k])

print('\nContinuum case\n')
pprint(simplify(Cdash))

#---------------------------------------------------------#
#                          I_hat                          #
#---------------------------------------------------------#

print('\nCompute Ihat for consistent modulus\n')

n = MatrixSymbol('n', 3, 3)

Ihat = zeros(6, 6)
I_outer_I = zeros(6, 6)
I_dev = zeros(6, 6)
n_outer_n = zeros(6, 6)
I4 = zeros(6, 6)

C_e = zeros(6, 6)

λ, gamma, beta, mu = symbols('λ, gamma, beta, mu')

for i in range(0, 6):
    a, b = voigt_map[i]
    for j in range(0, 6):
        c, d = voigt_map[j]

        δab = S(int(a == b))
        δac = S(int(a == c))
        δad = S(int(a == d))
        δbc = S(int(b == c))
        δcd = S(int(c == d))
        δbd = S(int(b == d))

        I_outer_I[i, j] = δab * δcd
        n_outer_n[i, j] = n[a, b] * n[c, d]

        I4[i, j] = Rational(1, 2) * ( δac*δbd + δad*δbc )

        I_dev[i, j] = I4[i, j] - Rational(1, 3) * I_outer_I[i, j]

        # Continuum mechanics tensors

        Ihat[i, j] = I_dev[i, j] - n_outer_n[i, j]

        C_e[i, j] = λ * I_outer_I[i, j] + 2 * mu * I4[i, j]

pprint(Ihat)

#--------------------------------------------------------------#
#               Algorithmically consistent modulus             #
#--------------------------------------------------------------#


Calg = (λ +  Rational(2, 3) * mu) * I_outer_I + 2 * mu * beta * I_dev # - 2 * mu * (gamma - 1 + beta) * n_outer_n

pprint(factor(Calg))
