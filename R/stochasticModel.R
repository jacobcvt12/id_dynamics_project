stochastic.dx.dt <- function(step.size, y, param) {
    ## probabilities

    # antibiotic treatment
    p.anti.tx <- 1 - exp(-step.size * param["alpha"] * y["R"])

    # restoration of colonization resistance
    p.col.res <- 1 - exp(-step.size * param["theta"] * y["S"])

    # treatment success
    p.tx.suc <- 1 - exp(-step.size * param["p"] * param["epsilon"])

    # col w/o immune response
    p.col.minus <- 1 - exp(-step.size * (1 - param["f"]) * param["lambda"] *
                           y["S"])

    # col w/ immune reponse
    p.col.plus <- 1 - exp(-step.size * param["f"] * param["lambda"] * y["S"])

    # disease
    p.disease <- 1 - exp(-step.size * param["phi"] * y["C"])

    # discharge rates
    p.res.dis <- 1 - exp(-step.size * 1 / param["k.r"])
    p.sus.dis <- 1 - exp(-step.size * 1 / param["k"])
    p.c.min.dis <- p.c.plu.dis <- p.sus.dis
    p.dis.dis <- 1 - exp(-step.size * 1 / param["k.d"])

    # calculate changes across compartments
    anti.tx <- rbinom(1, y["R"], p.anti.tx)
    col.rs <- rbinom(1, y["S"], p.col.res)
    tx.suc <- rbinom(1, y["D"], p.tx.suc)
    col.minus <- rbinom(1, y["S"], p.col.minus)
    col.plus <- rbinom(1, y["S"], p.col.plus)
    disease <- rbinom(1, y["C"], p.disease

    res.dis <- rbinom(1, y["R"], p.res.dis)
    sus.dis <- rbinom(1, y["S"], p.sus.dis)
    c.min.dis <- rbinom(1, y["C"], p.c.min.dis)
    c.plu.dis <- rbinom(1, y["P"], p.c.plu.dis)
    dis.dis <- rbinom(1, y["D"], p.dis.dis)

    # deltas for compartments
    dR <- -anti.tx + col.rs - res.dis
    dS <- anti.tx - col.rs - sus.dis
    dC <- col.minus - disease - c.min.dis
    dP <- col.plus - c.plu.dis
    dD <- disease - tx.suc - dis.dis

    return(c(dR, dS, dC, dP, dD))
}
