stochastic.sim <- function(param, initial.state,
                           dx.dt.func,
                           max.step=1e3,
                           iters=1e3) {

    # initialize matrix for output
    sim.data <- matrix(ncol=2, nrow=iters)

    # run simulation
    for (s in 1:iters) {
        sim.data[s, ] <- stochastic.run(param, initial.state,
                                        dx.dt.func, max.step)
    }

    # calculate 95% quantile based interval for attack rate
    attack.rate <- sim.data[, 2] / sim.data[, 1]
    sim.ci <- quantile(attack.rate, c(0.025, 0.5, 0.975))

    return(sim.ci)
}

stochastic.run <- function(param, initial.state, 
                           dx.dt.func, 
                           max.step) {
    # let's assume the model is NOT freq dependent for simplicity

    # initialize time, y, and counters for pts and diseased
    t <- 1
    y <- initial.state

    # consider are "susceptibles" to be the resistant and 
    # to start
    uninf <- sum(y[c("R", "S")])

    # any initial diseased people are not initially counted
    inf <- 0

    # store values
#     sir.output <- matrix(ncol=1 + length(y), nrow=max.step)
#     colnames(sir.output) <- c("time", names(y))
#     sir.output[1,] <- c(t, y)

    # while there are no colonized (w/ or w/o protection) or diseased
    # and while we are lt max step
    while (sum(y[c("C", "P", "D")], na.rm=TRUE) & t < max.step) { 
        # NAs removed for models that collapse colonized to one compartment
        t <- t + 1

        # dx.dt.func should be passed as function name 
        # not as a string
        delta <- dx.dt.func(1, y, param)

        y <- y + delta$changes
        inf <- inf + delta$colonizations
        uninf <- uninf + delta$resistant.again

#         sir.output[t, ] <- c(t, y)
    }

#     sir.output <- sir.output[rowSums(is.na(sir.output)) != 
#                              ncol(sir.output), ]
#     sir.output <- na.omit(sir.output)
    
#     return(sir.output)
    return(c(uninfected=uninf, infected=inf))
}

stochastic.dx.dt <- function(step.size, y, param) {
    R <- y["R"]
    S <- y["S"]
    P <- y["P"]
    C <- y["C"]
    D <- y["D"]

    # calculate lambda
    lambda <- param["beta.c"] * (C + P) + 
              param["beta.d"] * D
    ## probabilities

    # antibiotic treatment
    p.anti.tx <- 1 - exp(-step.size * param["alpha"])

    # restoration of colonization resistance
    p.col.res <- 1 - exp(-step.size * param["theta"])

    # treatment success
    p.tx.suc <- 1 - exp(-step.size * param["p"] * param["epsilon"])

    # col w/o immune response
    p.col.minus <- 1 - exp(-step.size * (1 - param["f"]) * lambda * S)

    # col w/ immune reponse
    p.col.plus <- 1 - exp(-step.size * param["f"] * lambda * S)

    # disease
    p.disease <- 1 - exp(-step.size * param["phi"] * C)

    # discharge rates
    p.res.dis <- 1 - exp(-step.size * param["k.r"])
    p.sus.dis <- 1 - exp(-step.size * param["k"])
    p.c.min.dis <- p.c.plu.dis <- p.sus.dis
    p.dis.dis <- 1 - exp(-step.size * param["k.d"])

    # calculate changes across compartments
    anti.tx <- rbinom(1, R, p.anti.tx)
    col.rs <- rbinom(1, S, p.col.res)
    tx.suc <- rbinom(1, D, p.tx.suc)
    col.minus <- rbinom(1, S, p.col.minus)
    col.plus <- rbinom(1, S, p.col.plus)
    disease <- rbinom(1, C, p.disease)

    res.dis <- rbinom(1, R, p.res.dis)
    sus.dis <- rbinom(1, S, p.sus.dis)
    c.min.dis <- rbinom(1, C, p.c.min.dis)
    c.plu.dis <- rbinom(1, P, p.c.plu.dis)
    dis.dis <- rbinom(1, D, p.dis.dis)

    # deltas for compartments
    dR <- -anti.tx + col.rs - res.dis
    dS <- anti.tx - col.rs - sus.dis - col.minus - col.plus + tx.suc
    dC <- col.minus - disease - c.min.dis
    dP <- col.plus - c.plu.dis
    dD <- disease - tx.suc - dis.dis

    # make sure none of the deltas are bigger than the prev comp size
    dR <- max(dR, -R)
    dS <- max(dS, -S)
    dC <- max(dC, -C)
    dP <- max(dP, -P)
    dD <- max(dD, -D)

    # careful! need to return in same order as y
    return(list(changes=c(dR, dS, dP, dC, dD),
                colonizations=(col.minus + col.plus),
                resistant.again=col.rs))
}

    stochastic.dx.dt.new <- function(step.size, y, param) {
      R <- y["R"]
      S.abx <- y["S.abx"]
      S.ft <- y["S.ft"]
      C <- y["C"]
      D <- y["D"]
      
      # calculate lambda
      lambda <- param["beta.c"] * C + 
        param["beta.d"] * D
      ## probabilities
      
      # antibiotic treatment in resistant
      p.anti.tx <- 1 - exp(-step.size * param["alpha"])
      
      # restoration of colonization resistance after antibiotics
      p.col.res.abx <- 1 - exp(-step.size * param["theta.abx"])
      
      # restoration of colonization resistance after fecal transplant
      p.col.res.ft <- 1 - exp(-step.size * param["theta.ft"])
      
      # antibiotic treatment success (in diseased)
      p.tx.suc.abx <- 1 - exp(-step.size * param["p"] * param["epsilon"] * (1 - param["tau"]))
      
      # fecal transplant treatment success
      p.tx.suc.ft <- 1 - exp(-step.size * param["rho"] * param["e"] * param["tau"])
      
      # col after antibiotics
      p.col.abx <- 1 - exp(-step.size * lambda * S.abx)
      
      # col after fecal transplant
      p.col.ft <- 1 - exp(-step.size * lambda * S.ft)
      
      # disease in colonized
      p.disease <- 1 - exp(-step.size * param["phi"] * C)
      
      # discharge rates
      p.res.dis <- 1 - exp(-step.size * param["k.r"])
      p.sus.ft.dis <- 1 - exp(-step.size * param["k"])
      p.col.dis <- p.sus.abx.dis <- p.sus.ft.dis
      p.dis.dis <- 1 - exp(-step.size * param["k.d"])
      
      # calculate changes across compartments
      anti.tx <- rbinom(1, R, p.anti.tx)
      col.res.abx <- rbinom(1, S.abx, p.col.res.abx)
      col.res.ft <- rbinom(1, S.ft, p.col.res.ft)
      tx.suc.abx <- rbinom(1, D, p.tx.suc.abx)
      tx.suc.ft <- rbinom(1, D, p.tx.suc.ft)
      col.abx <- rbinom(1, S.abx, p.col.abx)
      col.ft <- rbinom(1, S.ft, p.col.ft)
      disease <- rbinom(1, C, p.disease)
      
      res.dis <- rbinom(1, R, p.res.dis)
      sus.abx.dis <- rbinom(1, S, p.sus.abx.dis)
      sus.ft.dis <- rbinom(1, S, p.sus.ft.dis)
      col.dis <- rbinom(1, C, p.col.dis)
      dis.dis <- rbinom(1, D, p.dis.dis)
      
      # deltas for compartments
      dR <- -anti.tx + col.res.abx + col.res.ft - res.dis
      dS.abx <- anti.tx - col.res.abx - sus.abx.dis - col.abx + tx.suc.abx
      dS.ft <- col.res.ft - sus.ft.dis - col.ft + tx.suc.ft
      dC <- col.abx + col.ft - disease - col.dis
      dD <- disease - tx.suc.abx - tx.suc.ft - dis.dis
      
      # make sure none of the deltas are bigger than the prev comp size
      dR <- max(dR, -R)
      dS.abx <- max(dS.abx, -S.abx)
      dS.ft <- max(dS.ft, -S.ft)
      dC <- max(dC, -C)
      dD <- max(dD, -D)
      
      # careful! need to return in same order as y
      return(list(changes=c(dR, dS.abx, dS.ft, dC, dD),
                  colonizations=(col.abx + col.ft),
                  resistant.again=col.res.abx + col.res.ft))
}
    
    stochastic.dx.dt.new <- function(step.size, y, param) {
      R <- y["R"]
      S.abx <- y["S.abx"]
      S.ft <- y["S.ft"]
      C <- y["C"]
      D <- y["D"]
      
      # calculate lambda
      lambda <- param["beta.c"] * C + 
        param["beta.d"] * D
      ## probabilities
      
      # antibiotic treatment in resistant
      p.anti.tx <- 1 - exp(-step.size * param["alpha"])
      
      # restoration of colonization resistance after antibiotics
      p.col.res.abx <- 1 - exp(-step.size * param["theta.abx"])
      
      # restoration of colonization resistance after fecal transplant
      p.col.res.ft <- 1 - exp(-step.size * param["theta.ft"])
      
      # antibiotic treatment success (in diseased)
      p.tx.suc.abx <- 1 - exp(-step.size * param["p"] * param["epsilon"] * (1 - param["tau"]))
      
      # fecal transplant treatment success
      p.tx.suc.ft <- 1 - exp(-step.size * param["rho"] * param["e"] * param["tau"])
      
      # col after antibiotics
      p.col.abx <- 1 - exp(-step.size * lambda * S.abx)
      
      # col after fecal transplant
      p.col.ft <- 1 - exp(-step.size * lambda * S.ft)
      
      # disease in colonized
      p.disease <- 1 - exp(-step.size * param["phi"] * C)
      
      # discharge rates
      p.res.dis <- 1 - exp(-step.size * param["k.r"])
      p.sus.ft.dis <- 1 - exp(-step.size * param["k"])
      p.col.dis <- p.sus.abx.dis <- p.sus.ft.dis
      p.dis.dis <- 1 - exp(-step.size * param["k.d"])
      
      # calculate changes across compartments
      anti.tx <- rbinom(1, R, p.anti.tx)
      col.res.abx <- rbinom(1, S.abx, p.col.res.abx)
      col.res.ft <- rbinom(1, S.ft, p.col.res.ft)
      tx.suc.abx <- rbinom(1, D, p.tx.suc.abx)
      tx.suc.ft <- rbinom(1, D, p.tx.suc.ft)
      col.abx <- rbinom(1, S.abx, p.col.abx)
      col.ft <- rbinom(1, S.ft, p.col.ft)
      disease <- rbinom(1, C, p.disease)
      
      res.dis <- rbinom(1, R, p.res.dis)
      sus.abx.dis <- rbinom(1, S, p.sus.abx.dis)
      sus.ft.dis <- rbinom(1, S, p.sus.ft.dis)
      col.dis <- rbinom(1, C, p.col.dis)
      dis.dis <- rbinom(1, D, p.dis.dis)
      
      # deltas for compartments
      dR <- -anti.tx + col.res.abx + col.res.ft - res.dis
      dS.abx <- anti.tx - col.res.abx - sus.abx.dis - col.abx + tx.suc.abx
      dS.ft <- col.res.ft - sus.ft.dis - col.ft + tx.suc.ft
      dC <- col.abx + col.ft - disease - col.dis
      dD <- disease - tx.suc.abx - tx.suc.ft - dis.dis
      
      # make sure none of the deltas are bigger than the prev comp size
      dR <- max(dR, -R)
      dS.abx <- max(dS.abx, -S.abx)
      dS.ft <- max(dS.ft, -S.ft)
      dC <- max(dC, -C)
      dD <- max(dD, -D)
      
      # careful! need to return in same order as y
      return(list(changes=c(dR, dS.abx, dS.ft, dC, dD),
                  colonizations=(col.abx + col.ft),
                  resistant.again=col.res.abx + col.res.ft))
}
