#! /usr/bin/env python
from sympy import diff, latex, symbols
from sympy.matrices import Matrix

# symbols
# parameters
a_c, a_d, beta_c, beta_d, delta, e, epsilon = \
    symbols('a_c, a_d, beta_c, beta_d, delta, e, epsilon')
k, k_d, lambda_sym, p, psi, rho, tau = \
    symbols('k, k_d, lambda_sym, p, psi, rho, tau')

# compartments
C, S_a, S_f, D, N = symbols('C, S_a, S_f, D, N')

# create elements of the "Fancy" F and V matrices as
# identifed in ODEs
# rate of appearance of new infections in compartment
F_fancy_1 = beta_c * C * (S_a + S_f) + beta_d * D * (S_a + S_f)
F_fancy_2 = 0

# rate of transfer of individuals by other means (V^- - V^+)
V_fancy_1 = (psi * C - k * C) - (a_c * delta * N)
V_fancy_2 = (p * epsilon * D * (1 - tau) +
             rho * e * D * tau + k_d * D) - \
            (a_d * delta * N + psi * C)

# Jacobian of "Fancy" F and V
F = Matrix([[diff(F_fancy_1, C), diff(F_fancy_1, D)],
            [diff(F_fancy_2, C), diff(F_fancy_2, D)]])

V = Matrix([[diff(V_fancy_1, C), diff(V_fancy_1, D)],
            [diff(V_fancy_2, C), diff(V_fancy_2, D)]])

# Next Generation Matrix calculation
A = F * V.inv()
R_0 = A.eigenvals()  # look at eigenvalues to find dominant

# print out latex of results
print("F")
print(latex(F))

print("V")
print(latex(V))

print("NGM")
print(latex(A))

print("eigenvalues")
print(latex(list(R_0)[0]))
print(latex(list(R_0)[1]))
