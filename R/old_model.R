# libraries
library(deSolve)

# parameters
a.r <- 0.75
a.s <- 0.22
a.cn <- a.cp <- 0.01
a.d <- 0.01
alpha <- 0.5
theta <- 0.033
beta.c <- beta.d <- 0.007
f <- 0.6
epsilon <- 0.1
p <- 0.8
phi <- 0.2
k.r <- 0.33
k <- 0.15
k.d <- 0.068

# parameter vector
param <- c(a.r=a.r, a.s=a.s, a.cn=a.cn, a.cp=a.cp, a.d=a.d, alpha=alpha,
           theta=theta, beta.c=beta.c, beta.d=beta.d, f=f, 
           epsilon=epsilon, p=p, phi=phi, k.r=k.r, k=k, k.d=k.d)

# sequence of times
max.time <- 1e3
times <- seq(0, max.time, 1) 

# initial state
initial.state <- c("R"=20,
                   "S"=0,
                   "P"=0,
                   "C"=0,
                   "D"=0)

# diff eq function
dx.dt <- function(t, y, param) {
    # current state
    R <- y["R"]
    S <- y["S"]
    P <- y["P"]
    C <- y["C"]
    D <- y["D"]

    # params see lines 5 - 18
    a.r <- param["a.r"]
    a.s <- param["a.s"]
    a.cn <- a.cp <- param[""] ##not sure here?
    a.d <- param["a.d"]
    alpha <- param["alpha"]
    theta <- param["theta"]
    beta.c <- beta.d <- param[""] ##not sure here?
    f <- param["f"]
    epsilon <- param["epsilon"]
    p <- param["p"]
    phi <- param["phi"]
    k.r <- param["k.r"]
    k <- param["k"]
    k.d <- param["k.d"]
    
    # ODEs
    lambda <- param["beta.c"] * (C + P) + 
              param["beta.d"] * D
    N = R + S + C + P + D
    dR <- a.r * delta*N + theta*S - k.r*R - alpha * R
    dS <- a.s * delta * N + alpha * R + p * epsilon * D - 
          theta * S - k * S - lambda * S
    dP <- a.cp * delta * N + f * lambda * S - k * P
    dC <- a.cn * delta * N + (1 - f) * lambda * S - phi * C - k * C
    dD <- a.d * delta * N + phi * C - p * epsilon * D - k.d * D

    return(list(c(dR, dS, dP, dC, dD)))
}
