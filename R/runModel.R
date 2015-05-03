runModel <- function(a.r = 0.75,
         a.s = 0.22,
         a.cn = a.cp <- 0.01,
         a.d = 0.01,
         alpha = 0.5,
         theta = 0.033,
         beta.c = 0.007, 
         beta.d = 0.007,
         f = 0.6,
         epsilon = 0.1,
         e = 0.01,
         p = 0.8,
         rho = 0.8,
         phi = 0.2,
         k.r = 0.33,
         k = 0.15,
         k.d = 0.068,
         initial.state= c(R=1e3, S=0, C=1e1, D=0), #added pop2 and betas code here
         max.time = 25) {
  
  
  
  
  
  
  #create the parameter vector.
  param <- c(a.r=a.r, a.s=a.s, a.cn=a.cn, a.cp=a.cp, a.d=a.d, alpha=alpha,
             theta=theta, beta.c=beta.c, beta.d=beta.d, f=f, 
             epsilon=epsilon, e=e, p=p, rho=rho,  phi=phi, k.r=k.r, k=k, k.d=k.d) #added 4 betas and m code here
  
  #Sequence of times at which we want estimates..
  #   ..here we say daily until max.time
  times <- seq(0,max.time,1)
  
  #Run the ODE-Solver
  sir.output <- lsoda(initial.state, times, dx.dt.new, param) 
  
  #now we will return the output.
  return(sir.output)
  
  results.CDiff <- runModel()
  plot.ode(results.CDiff)
  legend.ode(15,800,results.CDiff)
}