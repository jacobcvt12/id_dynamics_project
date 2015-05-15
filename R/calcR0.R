calc.R0 <- function(N, param) {
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
  
  # calculations for R0 at disease-free equilibrium
  S.f.0 <- 0
  S.a.0 <- ((a.s.abx * k.r + alpha) * N) / (alpha + a.s.abx * k.r + (1 - a.s.abx) * k + theta.abx)
  R.0 <- -((S.a.0 + S.f.0) * (beta.c * epsilon.abx * p.abx * psi - 
                              beta.c * epsilon.abx * p.abx - 
                              beta.c * epsilon.ft * psi * p.ft + 
                              beta.c * k.d - 
                              beta.d * phi)) / 
         ((k - phi) * (epsilon.abx * p.abx * psi - 
                       epsilon.abx * p.abx - 
                       epsilon.ft * psi * p.ft - 
                       k.d))
  
  return(abs(R.0))
}
