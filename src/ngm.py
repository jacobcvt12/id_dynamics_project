#! /usr/bin/env python
from sympy import latex, symbols
from sympy.matrices import Matrix

# symbols
# parameters
a_c, a_d, beta_c, beta_d, delta, e, epsilon = \
    symbols('a_c, a_d, beta_c, beta_d, delta, e, epsilon')
k, k_d, lambda_sym, p, psi, rho, tau = \
    symbols('k, k_d, lambda_sym, p, psi, rho, tau')

# compartments
C, S_a, S_f, D, N = symbols('C, S_a, S_f, D, N')

# create elements of the "Fancy" (cursive) F and V matrices as
# identifed in ODEs
# rate of appearance of new infections in compartment
F_fancy = Matrix([beta_c * C * (S_a + S_f) + beta_d * D * (S_a + S_f), 0])

# rate of transfer of individuals by other means (V^- - V^+)
V_fancy = Matrix([(psi * C - k * C) - (a_c * delta * N),
                  (p * epsilon * D * (1 - tau) +
                   rho * e * D * tau + k_d * D) - (a_d * delta * N + psi * C)])

# Jacobian of "Fancy" F and V
F = F_fancy.jacobian([C, D])
V = V_fancy.jacobian([C, D])

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
