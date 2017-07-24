
import numpy as np

from sympy import *
from sympy import eye
from sympy.tensor.array import Array, MutableDenseNDimArray, tensorproduct, tensorcontraction
from sympy.tensor.tensor import *

init_printing()

voigt_map = [[0, 0], [1, 1], [2, 2], [1, 2], [0, 2], [0, 1]]

voigt_dict = {(0, 0) : 0,
              (1, 1) : 1,
              (2, 2) : 2,
              (1, 2) : 3,
              (2, 1) : 3,
              (2, 0) : 4,
              (0, 2) : 4,
              (0, 1) : 5,
              (1, 0) : 5}

def IndexToVoigt(Cindex):

    C = MutableDenseNDimArray(range(6 * 6), (6, 6))

    for a in range(0, 6):
        i, j = voigt_map[a]
        for b in range(0, 6):
            k, l = voigt_map[b]
            C[a, b] = Cindex[i, j, k, l]
    return C

def Identity4():
    FourthI = MutableDenseNDimArray(range(81), (3, 3, 3, 3))
    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            for k in [0, 1, 2]:
                for l in [0, 1, 2]:
                    δik = int(i == k)
                    δjl = int(j == l)
                    δil = int(i == l)
                    δjk = int(j == k)

                    FourthI[i,j,k,l] = Rational(1, 2) * (δik * δjl + δil * δjk)

    return FourthI

def DeviatoricProjection():
    """ P = I - 1/3 * 1⊗1 """
    devProj = MutableDenseNDimArray(range(3*3*3*3), (3, 3, 3, 3))
    I = Identity4()
    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            for k in [0, 1, 2]:
                for l in [0, 1, 2]:
                    δij = int(i == j)
                    δik = int(i == k)
                    δil = int(i == l)
                    δjk = int(j == k)
                    δkl = int(k == l)
                    δjl = int(j == l)

                    devProj[i,j,k,l] = I[i,j,k,l] - Rational(1, 3) * δij*δkl
    return devProj

def P_ddot_D_ddot_P(P, D):
    D_dash = MutableDenseNDimArray( zeros(81), (3,3,3,3) )

    for i in range(3):
         for j in range(3):
                 for k in range(3):
                         for l in range(3):
                                 for m in range(3):
                                         for n in range(3):
                                                 for o in range(3):
                                                         for p in range(3):
                                                                 D_dash[i,j,o,p] += P[i,j,k,l]*D[k,l,m,n]*P[m,n,o,p]
    return D_dash

def SymbolicCbar():
    C_bar = symarray('', (3,3,3,3) )
    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            for k in [0, 1, 2]:
                for l in [0, 1, 2]:
                    C_bar[i,j,k,l] = symbols('C_bar(' + str(voigt_dict[(i, j)]) + '\,' + str(voigt_dict[(k, l)]) + ')')
    C_bar = Array(C_bar, (3, 3, 3, 3))
    return C_bar

I = Identity4()

τ00 = Symbol('τ(0, 0)')
τ11 = Symbol('τ(1, 1)')
τ22 = Symbol('τ(2, 2)')
τ12 = Symbol('τ(1, 2)')
τ02 = Symbol('τ(0, 2)')
τ01 = Symbol('τ(0, 1)')
τ = Matrix(([τ00, τ01, τ02], [τ01, τ11, τ12], [τ02, τ12, τ22]))

C_bar = SymbolicCbar()

# The D in P : D : P
D = C_bar + Rational(2, 3) * tensorcontraction( τ, (0, 1) ) * I \
          - Rational(2, 3) * ( tensorproduct( τ, eye(3) ) + tensorproduct( eye(3), τ ) )

D_dash = P_ddot_D_ddot_P(DeviatoricProjection(), D)

κ, p = symbols('κ, p')

C = (κ + p) * tensorproduct(eye(3), eye(3)) - 2 * p * I + D_dash

D_dash_voigt = IndexToVoigt(D_dash)

for i in range(6):
    for j in range(6):
        D_dash_voigt[i, j] = simplify(D_dash_voigt[i, j])

print('Material tangent matrix in Voigt notation')
print(D_dash_voigt)
