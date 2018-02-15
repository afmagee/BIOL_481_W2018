sampleLKtrajectory <- function(H0,P0,alpha1,alpha2,beta1,beta2,sampled.times,n.steps) {
  # recover()
  
  # Calculate time slice size
  h <- (max(sampled.times) - min(sampled.times))/(n.steps)
  
  # For storing information, nsteps+1 because we're starting at t=0
  Ht <- numeric(n.steps+1)
  Pt <- numeric(n.steps+1)
  
  # Given initial sizes
  Ht[1] <- H0
  Pt[1] <- P0
  
  # A step at a time, simulate next step
  for (t in 1:n.steps) {
    Ht[t+1] <- Ht[t]+(alpha1 - beta1*Pt[t])*Ht[t]*h
    Pt[t+1] <- Pt[t]+(-alpha2 + beta2*Ht[t])*Pt[t]*h
  }
  
  # Filter to only the ones we needed
  keep <- (min(sampled.times)+(0:n.steps)*h) %in% sampled.times
  
  if ( sum(keep) != length(sampled.times) ) {
    stop("Error in sampleLKtrajectory, simulated times not matching real times. Are times evenly spaced?")
  }
  
  Ht <- Ht[keep]
  Pt <- Pt[keep]
  
  return(data.frame(Ht=Ht,Pt=Pt))
  
}

# reparamLN <- function(mean,sd) {
#   meanlog <- log(mean/sqrt(1+sd/(mean^2)))
#   sdlog <- sqrt(log(1+sd/(mean^2)))
#   return(data.frame(meanlog=meanlog,sdlog=sdlog))
# }

logLikelihoodLK <- function(data,
                            log_H0,
                            log_P0,
                            alpha1,
                            alpha2,
                            beta1,
                            beta2,
                            log_sigma1,
                            log_sigma2,
                            n.steps) {
  # recover()
  
  # un-log parameters
  # We have them input on the log scale so optimization is unbounded
  H0 <- exp(log_H0)
  P0 <- exp(log_P0)
  sigma1 <- exp(log_sigma1)
  sigma2 <- exp(log_sigma2)
  
  # Need to know what we think the pop sizes were at the times we sampled
  path <- sampleLKtrajectory(H0=H0,P0=P0,alpha1=alpha1,alpha2=alpha2,beta1=beta1,beta2=beta2,sampled.times=data$time,n.steps=n.steps)
  
  # Now it's just a measurement model
  lnL <- sum(dlnorm(data$prey,log(path$Ht),sigma1,log=TRUE)) + # measurements of prey pop sizes
    sum(dlnorm(data$pred,log(path$Pt),sigma2,log=TRUE)) # measurements of pred pop sizes
  
  return(lnL)
  
}
                          
fitLotkaVolteraMaximumLikelihood <- function(data,...) {
  recover()
  
  # Make sure we have recognizable names in data
  if ( !any(grepl("time",names(data))) ) {
    stop("data must be data frame with a column named time or times")
  }
  if ( !any(grepl("pred",names(data))) ) {
    stop("data must be data frame with a column named predator or unambiguous abbreviation thereof")
  }
  if ( !any(grepl("prey",names(data))) ) {
    stop("data must be data frame with a column named prey")
  }
  
  ## First-pass estimation with derivatives, following Howard
  
  # dH/dt = (alpha1 - beta1*P)*H
  # 1/H * dH/dt = alpha1 - beta1*P
  # dP/dt = (-alpha2 + beta2*H)*P
  # 1/P * dP/dt = alpha2 - beta2*H
  # We know everything but the alphas/betas so we estimate them by least squares (regression)
  # This is a passable estimate, but dP/dt is already an estimate, so there's extra error
  # Instead, we start here and optimize with likelihood
  dH <- (data$prey[3:dim(data)[1]] - data$prey[1:(dim(data)[1]-2)])/(data$time[3:dim(data)[1]] - data$time[1:(dim(data)[1]-2)]) #nrow is slower than dim()[1] for esoteric reasons
  dP <- (data$pred[3:dim(data)[1]] - data$pred[1:(dim(data)[1]-2)])/(data$time[3:dim(data)[1]] - data$time[1:(dim(data)[1]-2)])
  
  dH_over_H <- dH/data$prey[2:(dim(data)[1]-1)]
  dP_over_P <- dP/data$pred[2:(dim(data)[1]-1)]
  
  # fit the lines
  least_squares_1 <- lm(dH_over_H ~ data$pred[2:(dim(data)[1]-1)])
  least_squares_2 <- lm(dP_over_P ~ data$prey[2:(dim(data)[1]-1)])
  
  alpha1_start <- as.numeric(least_squares_1$coefficients[1]) # we wrap these with as.numeric to get rid of obnoxious names
  alpha2_start <- as.numeric(-least_squares_2$coefficients[1]) # The regression estimate assumes we are adding all the parameters, but this is subracted, so we take the negative
  
  beta1_start <- as.numeric(-least_squares_1$coefficients[2]) # The regression estimate assumes we are adding all the parameters, but this is subracted, so we take the negative
  beta2_start <- as.numeric(least_squares_2$coefficients[2])
  
  # recover()
  
  # Pick a time step size
  # This default will work reasonably well for data with two peaks
  # We want about 4000-5000 steps
  total_time <- max(data$time) - min(data$time)
  
  # For final optimization, smaller steps
  steps_per_interval_final <- floor(5000/total_time)
  n_steps_final <- total_time*steps_per_interval_final
  
  # For initial, less precise is fine
  steps_per_interval_initial <- floor(1000/total_time)
  n_steps_initial <- total_time*steps_per_interval_initial
  
  # Find area near peak for starting sizes and errors, conditioning on approximate Lotka-Voltera parameters
  # If we don't do this, we can easily get lost in crazy parameters (optimization is unstable!)
  # A imore principled approach is to start optimization many times in many places (parameter combinations), but this is too slow for a 2 hour lab, given the time to optimize
  fnLKInitAndVar <- function(par) {
    H0 <- par[1]
    P0 <- par[2]
    alpha1 <- alpha1_start
    alpha2 <- alpha2_start
    beta1 <- beta1_start
    beta2 <- beta2_start
    sigma1 <- par[3]
    sigma2 <- par[4]
    -logLikelihoodLK(data,H0,P0,alpha1,alpha2,beta1,beta2,sigma1,sigma2,n.steps=n_steps_initial)
  }
  
  opts <- optim(c(log(data$pred[1]),log(data$prey[1]),log(0.3),log(0.3)),fnLKInitAndVar,method="BFGS")
  
  fnLKforOptimML <- function(par) {
    H0 <- par[1]
    P0 <- par[2]
    alpha1 <- par[3]
    alpha2 <- par[4]
    beta1 <- par[5]
    beta2 <- par[6]
    sigma1 <- par[7]
    sigma2 <- par[8]
    -logLikelihoodLK(data,H0,P0,alpha1,alpha2,beta1,beta2,sigma1,sigma2,n.steps=n_steps_final)
  }
  opt_ml <- optim(c(opts$par[1:2],alpha1_start,alpha2_start,beta1_start,beta2_start,opts$par[3:4]),fnLKforOptimML,method="BFGS",hessian=TRUE)
  
  ml_par <- opt_ml$par
  ml_par[c(1:2,7:8)] <- exp(ml_par[c(1:2,7:8)]) # move off log-scale
  names(ml_par) <- c("H0","P0","alpha1","alpha2","beta1","beta2","sigma1","sigma2")
  
  # recover()
  
  # Compute confidence intervals
  asymptotic_sd <- sqrt(1/diag(opt_ml$hessian))
  asymptotic_025 <- numeric(8)
  asymptotic_025[c(1:2,7:8)] <- qlnorm(0.025,opt_ml$par[c(1:2,7:8)],asymptotic_sd[c(1:2,7:8)])
  asymptotic_025[c(3:6)] <- qnorm(0.025,opt_ml$par[c(3:6)],asymptotic_sd[3:6])
  asymptotic_975 <- numeric(8)
  asymptotic_975[c(1:2,7:8)] <- qlnorm(0.975,opt_ml$par[c(1:2,7:8)],asymptotic_sd[c(1:2,7:8)])
  asymptotic_975[c(3:6)] <- qnorm(0.975,opt_ml$par[c(3:6)],asymptotic_sd[3:6])
  
  asymptotic_ci <- cbind(asymptotic_025,asymptotic_975)
  colnames(asymptotic_ci) <- c("2.5%","97.5%")
  rownames(asymptotic_ci) <- c("H0","P0","alpha1","alpha2","beta1","beta2","sigma1","sigma2")
  
  return(list(logLikelihood=-opt_ml$value,maximum.likelihood.parameter.estimates=ml_par,confidence.intervals=asymptotic_ci))
  
}

assessLKModelFit <- function(data,ml.parameters,n.sim=100,...) {
  if ( hasArg(pred.color) ) {
    pred_color <- list(...)$pred.color
  } else {
    pred_color <- "#cc000090"
  }
  if ( hasArg(prey.color) ) {
    prey_color <- list(...)$prey.color
  } else {
    prey_color <- "#3366ff90"
  }
  
  # recover()
  
  # Useful things
  n_time_points <- length(data$times)
  
  # Simulate trajectories from model
  pred_sim <- matrix(NA,nrow=n_time_points,ncol=n.sim)
  prey_sim <- matrix(NA,nrow=n_time_points,ncol=n.sim)
  pb <- txtProgressBar(min = 0, max = n.sim, style = 3) # report progress
  for (i in 1:n.sim) {
    # Simulate true trajectory
    sim_traj <- sampleLKtrajectory(H0=ml.parameters[1],P0=ml.parameters[2],alpha1=ml.parameters[3],alpha2=ml.parameters[4],beta1=ml.parameters[5],beta2=ml.parameters[6],t.sim=max(data$times),sampled.times=data$times,h=1e-3)
    # add measurement error
    prey_sim[,i] <- rlnorm(n_time_points,log(sim_traj$Ht),ml.parameters[7])
    pred_sim[,i] <- rlnorm(n_time_points,log(sim_traj$Pt),ml.parameters[8])
    setTxtProgressBar(pb, i) # report progress
  }
  
  # Plot model

  # # so we can put axes on same scale
  # plot_max <- max(c(max(prey_sim),max(pred_sim))) # since there's measurement error, it's possible that there's a super-high predator pop size
  
  par(mfrow=c(2,1),mai=c(0.8,0.8,0.2,0.2),xpd=TRUE)
  matplot(x=data$times,prey_sim,col=prey_color,type="l",lty=1,xlab="time",ylab="prey pop size")
  lines(data$times,data$prey,lwd=3)

  matplot(x=data$times,pred_sim,col=pred_color,type="l",lty=1,xlab="time",ylab="predator pop size")
  lines(data$times,data$pred,lwd=3)
  par(mfrow=c(1,1))
  
}

AGD <- function(fn,init,)