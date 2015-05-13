# delta * N is the number of new admits per step
# diff eq function
dx.dt.orig <- function(t, y, param, delta=0) {
    # current state
    R <- y["R"]
    S <- y["S"]
    P <- y["P"]
    C <- y["C"]
    D <- y["D"]

    # params 
    a.r <- param["a.r"]
    a.s <- param["a.s"]
    a.cn <- param["a.cn"]
    a.cp <- param["a.cp"]
    a.d <- param["a.d"]
    alpha <- param["alpha"]
    theta <- param["theta"]
    beta.c <- param["beta.c"]
    beta.d <- param["beta.d"] 
    f <- param["f"]
    epsilon <- param["epsilon"]
    p <- param["p"]
    phi <- param["phi"]
    k.r <- param["k.r"]
    k <- param["k"]
    k.d <- param["k.d"]
    
    lambda <- param["beta.c"] * (C + P) + 
              param["beta.d"] * D
    N = R + S + C + P + D

    # ODEs
    dR <- a.r * delta*N + theta*S - k.r*R - alpha * R
    dS <- a.s * delta * N + alpha * R + p * epsilon * D - 
          theta * S - k * S - lambda * S
    dP <- a.cp * delta * N + f * lambda * S - k * P
    dC <- a.cn * delta * N + (1 - f) * lambda * S - phi * C - k * C
    dD <- a.d * delta * N + phi * C - p * epsilon * D - k.d * D

    return(list(c(dR, dS, dP, dC, dD)))
}

dx.dt.new <- function(t, y, param, delta = 1e-1) {
    # current state
    R <- y["R"]
    S.a <- y["S.abx"]
    S.f <- y["S.ft"]
    C <- y["C"]
    D <- y["D"]

    # params 
    a.r <- param["a.r"]   #all a's are respective admission proportion (delta is overall admit rate/N/day)
    a.s.abx <- param["a.s.abx"]
    a.s.ft <- param["a.s.ft"]
    a.c <- param["a.c"]
    a.d <- param["a.d"]
    alpha <- param["alpha"]   # abx Rx rate in R (community) (R -> S rate)
    theta.abx <- param["theta.abx"]   #S.abx -> R rate
    theta.ft <- param["theta.ft"]   #S.ft -> R rate
    beta.c <- param["beta.c"]
    beta.d <- param["beta.d"]   
    epsilon <- param["epsilon"]   #abx Tx rate /day (for diseased)
    e <- param["e"]   #fecal trsplt rate /day (diseased)
    p <- param["p"]   #prob. of abx Tx success
    rho <- param["rho"]   #prob. of fecal transplant success
    tau <- param["tau"]   # proportion of patients (D) receiving fecal transplant 
    phi <- param["phi"]   # disease rate/day in C
    k.r <- param["k.r"]   # R dc rate /day
    k <- param["k"]   # S & C dc rate 
    k.d <- param["k.d"]   # D dc rate
    
    lambda <- param["beta.c"] * C + 
              param["beta.d"] * D    # transmission (S -> C)
    N = R + S.abx + S.ft + C + D

    # ODEs
    dR <- a.r * delta*N + theta.abx*S.abx + theta.ft * S.ft - k.r * R - alpha * R
    dS.ft <- rho * e * D * tau - 
           theta.ft * S.ft - k * S.ft - lambda * S.ft
    dS.abx <- a.s.abx * delta * N + alpha * R + p * epsilon * D * (1-tau) - 
           theta.abx * S.abx - k * S.abx - lambda * S.abx
    dC <- a.c * delta * N + lambda * S.abx + lambda * S.ft - phi * C - k * C
    dD <- a.d * delta * N + phi * C - p * epsilon * D * (1-tau) - rho * e * D * tau - k.d * D

    return(list(c(dR, dS, dC, dD)))
}
