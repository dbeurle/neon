
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
        
        δik = float(i == k)
        δjl = float(j == l)
        δil = float(i == l)
        δjk = float(j == k)
        
        Cdash[a, b] = 0.5 * (δik*σ[j,l] + δil*σ[j,k] + δjk*σ[i,l] + δjl*σ[i,k])

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
        
        δik = float(i == k)
        δjl = float(j == l)
        δil = float(i == l)
        δjk = float(j == k)
        
        Cdash[a, b] = 0.5 * (δik*σ[j,l] + δil*σ[j,k] + δjk*σ[i,l] + δjl*σ[i,k])

print('\nContinuum case\n')
pprint(Cdash)

#---------------------------------------------------------#
#                          I_hat                          #
#---------------------------------------------------------#

print('\nCompute Ihat for consistent modulus\n')

n = MatrixSymbol('n', 3, 3)

Ihat = zeros(6, 6)

for i in range(0, 6):
    a, b = voigt_map[i]
    for j in range(0, 6):
        c, d = voigt_map[j]

        δab = float(a == b)
        δac = float(a == c)
        δad = float(a == d)
        δbc = float(b == c)
        δcd = float(c == d)
        δbd = float(b == d)
        
        # Ihat = I - 1/3 * 1-outer-1
        Ihat[i, j] = 0.5 * ( δac*δbd + δac*δbc ) - 1.0 / 3.0 * δab * δcd - n[a, b] * n[c, d]

pprint(Ihat)


