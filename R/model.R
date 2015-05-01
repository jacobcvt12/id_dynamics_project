# diff eq function
dx.dt.orig <- function(t, y, param) {
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

    # delta * N is the number of new admits per step
    delta <- 0 # closed system for now

    # ODEs
    dR <- a.r * delta*N + theta*S - k.r*R - alpha * R
    dS <- a.s * delta * N + alpha * R + p * epsilon * D - 
          theta * S - k * S - lambda * S
    dP <- a.cp * delta * N + f * lambda * S - k * P
    dC <- a.cn * delta * N + (1 - f) * lambda * S - phi * C - k * C
    dD <- a.d * delta * N + phi * C - p * epsilon * D - k.d * D

    return(list(c(dR, dS, dP, dC, dD)))
}

dx.dt.new <- function(t, y, param) {
    # current state
    R <- y["R"]
    S <- y["S"]
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
    N = R + S + C + D

    # delta * N is the number of new admits per step
    delta <- 0 # closed system for now

    # ODEs
    dR <- a.r * delta*N + theta*S - k.r*R - alpha * R
    dS <- a.s * delta * N + alpha * R + p * epsilon * D - 
          theta * S - k * S - lambda * S
    dC <- a.cn * delta * N + (1 - f) * lambda * S - phi * C - k * C
    dD <- a.d * delta * N + phi * C - p * epsilon * D - k.d * D

    return(list(c(dR, dS, dC, dD)))
}
