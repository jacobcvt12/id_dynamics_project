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
