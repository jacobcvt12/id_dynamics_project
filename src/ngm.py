#! /usr/bin/env python
from sympy import Symbol, diff, latex
from sympy.matrices import Matrix

# symbols
a_c = Symbol("a_c")
a_d = Symbol("a_d")
beta_c = Symbol("beta_c")
beta_d = Symbol("beta_d")
delta = Symbol("delta")
e = Symbol("e")
epsilon = Symbol("epsilon")
k = Symbol("k")
k_d = Symbol("k_d")
lambda_sym = Symbol("lambda")
p = Symbol("p")
psi = Symbol("phi")
rho = Symbol("rho")
tau = Symbol("psi")
C = Symbol("C")
S_a = Symbol("S_a")
S_f = Symbol("S_f")
D = Symbol("D")
N = Symbol("N")

# fancy matrices
F_fancy_1 = beta_c * C * (S_a + S_f) + beta_d * D * (S_a + S_f)
F_fancy_2 = 0

V_fancy_1 = (psi * C - k * C) - (a_c * delta * N)
V_fancy_2 = (p * epsilon * D * (1 - tau) +
             rho * e * D * tau + k_d * D) - \
            (a_d * delta * N + psi * C)

# F and V matrices
F = Matrix([[diff(F_fancy_1, C), diff(F_fancy_1, D)],
            [diff(F_fancy_2, C), diff(F_fancy_2, D)]])

V = Matrix([[diff(V_fancy_1, C), diff(V_fancy_1, D)],
            [diff(V_fancy_2, C), diff(V_fancy_2, D)]])

# Next Generation Matrix
A = F * V.inv()
R_0 = A.eigenvals()

print("F")
print(latex(F))

print("V")
print(latex(V))

print("NGM")
print(latex(A))

print("eigenvalues")
print(latex(R_0))

# DFE
# a_r, delta, N, theta, S, k_r, R, alpha = \
#     symbols('a_r, delta, N, theta, S, k_r, R, alpha')
# a_s, p, epsilon, D, k, beta_c, C_min, C_plu, beta_d = \
#     symbols('a_s, p, epsilon, D, k, beta_c, C_min, C_plu, beta_d')
#
# sys1 = Eq(a_r * delta * N + theta * S - k_r * R - alpha * R, 0)
# sys2 = Eq(a_s * delta * N + alpha * R- theta * S - k * S, 0)
