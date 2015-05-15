calc.R0 <- function(t, y, param, delta = 1e-1) {
  # current state
  R <- y["R"]
  S.ft <- y["S.ft"]
  S.abx <- y["S.abx"]
  C <- y["C"]
  D <- y["D"]
  
  # params 
  a.r <- param["a.r"]   #all a's are respective admission proportion (delta is overall admit rate/N/day)
  a.s.abx <- param["a.s.abx"]
  a.c <- param["a.c"]
  a.d <- param["a.d"]
  alpha <- param["alpha"]   # abx Rx rate in R (community) (R -> S rate)
  theta.abx <- param["theta.abx"]   #S.abx -> R rate
  theta.ft <- param["theta.ft"]   #S.ft -> R rate
  beta.c <- param["beta.c"]
  beta.d <- param["beta.d"]   
  epsilon.abx <- param["epsilon.abx"]   #abx Tx rate /day (for diseased)
  epsilon.ft <- param["epsilon.ft"]   #fecal trsplt rate /day (diseased)
  p.abx <- param["p.abx"]   #prob. of abx Tx success
  p.ft <- param["p.ft"]   #prob. of fecal transplant success
  psi <- param["psi"]   # proportion of patients (D) receiving fecal transplant 
  phi <- param["phi"]   # disease rate/day in C
  k.r <- param["k.r"]   # R dc rate /day
  k <- param["k"]   # S & C dc rate 
  k.d <- param["k.d"]   # D dc rate
  
  lambda <- param["beta.c"] * C + 
    param["beta.d"] * D    # transmission (S -> C)
  N = R + S.abx + S.ft + C + D
  
  # calculations for R0 at disease-free equilibrium
  S.f.0 <- 0
  S.a.0 <- ((a.s.abx * k.r + alpha) * N) / (alpha + a.s.abx * k.r + (1 - a.s.abx) * k + theta.abx)
  R.0 <- -((S.a.0 + S.f.0) * (-beta.c * epsilon.ft * psi * p.ft + beta.c * epsilon.abx * p.abx * psi - beta.c * epsilon.abx * p.abx - beta.c * k.d - beta.c * phi)) / 
    ((k - phi) * (-epsilon.ft * psi * p.ft + epsilon.abx * p.abx * psi - epsilon.abx * p.abx - k.d))
  
  return(R.0)
}