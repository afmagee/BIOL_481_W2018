# Here we simulate a trajectory for the predator and prey population sizes
# theta is a vector containing parameters
# LKmodel is a function that takes in theta, and the current predator and prey pop sizes, and returns a vector c(dPred/dt,dPrey/dt)
# In other words, this function can plot the trajectory for any (deterministic) model of predator and prey population cycling that depends only on predator and prey population sizes
sampleGeneralizedLKTrajectory <- function(prey.initial.size,pred.initial.size,theta,LKmodel,sampled.times,n.steps) {
  # recover()
  
  # Calculate time slice size
  # We don't solve any differential equations here
  # We just take tiny steps along the population trajectory given the slope at an instant
  # If the step size is sufficiently small, this is a passable approximation of more complicated techniques
  h <- (max(sampled.times) - min(sampled.times))/(n.steps)
  
  # For storing information, nsteps+1 because we're starting at t=0
  prey_trajectory <- numeric(n.steps+1)
  pred_trajectory <- numeric(n.steps+1)
  
  # Given initial sizes
  prey_trajectory[1] <- prey.initial.size
  pred_trajectory[1] <- pred.initial.size

  # A step at a time, simulate next step
  for (t in 1:n.steps) {
    trajectory_slope <- LKmodel(prey.size=prey_trajectory[t],pred.size=pred_trajectory[t],theta=theta)
    prey_trajectory[t+1] <- prey_trajectory[t]+trajectory_slope[1]*h
    pred_trajectory[t+1] <- pred_trajectory[t]+trajectory_slope[2]*h
  }
  
  # Filter to only the ones we needed
  keep <- (min(sampled.times)+(0:n.steps)*h) %in% sampled.times
  
  if ( sum(keep) != length(sampled.times) ) {
    stop("Error in sampleGeneralizedLKTrajectory, simulated times not matching real times. Are times evenly spaced?")
  }
  
  prey_trajectory <- prey_trajectory[keep]
  pred_trajectory <- pred_trajectory[keep]
  
  return(data.frame(prey_trajectory=prey_trajectory,pred_trajectory=pred_trajectory))
  
}

# Calculates the likelihood for a generalized class of Lotka-Voltera models
# Takes the following argumens:
#      data: the data, as a dataframe with columns named time(s), pred(ator), and prey
#      prey.initial.size: the population size for the prey at the first time point
#      pred.initial.size: the population size for the predator at the first time point
#      theta: the parameters for the model of population size change
#      data: a function that takes in theta, and the current predator and prey pop sizes, and returns a vector c(dPred/dt,dPrey/dt)
#      sigma1: the value of sigma for the (lognormal) error of prey measurements
#      sigma2: the value of sigma for the (lognormal) error of predator measurements
#      parameters.are.log.scale: are the parameter values given on the log scale? (They will be for optimization)
logLikelihoodGeneralizedLK <- function(data,
                                       prey.initial.size,
                                       pred.initial.size,
                                       theta,
                                       LKmodel,
                                       sigma1,
                                       sigma2,
                                       n.steps,
                                       parameters.are.log.scale=FALSE) {
  # recover()
  
  # un-log parameters
  # We may have them input on the log scale so optimization is unbounded
  if ( parameters.are.log.scale ) {
    prey.initial.size <- exp(prey.initial.size)
    pred.initial.size <- exp(pred.initial.size)
    theta <- exp(theta)
    sigma1 <- exp(sigma1)
    sigma2 <- exp(sigma2)
  }
  
  # Need to know what we think the pop sizes were at the times we sampled
  path <- sampleGeneralizedLKTrajectory(prey.initial.size=prey.initial.size,pred.initial.size=pred.initial.size,theta=theta,LKmodel=LKmodel,sampled.times=data$time,n.steps=n.steps)

  # Now it's just a measurement model
  lnL <- sum(dlnorm(data$prey,log(path$prey_trajectory),sigma1,log=TRUE)) + # measurements of prey pop sizes
    sum(dlnorm(data$pred,log(path$pred_trajectory),sigma2,log=TRUE)) # measurements of pred pop sizes
  
  return(lnL)
  
}
 
fitGeneralizedLotkaVolteraMaximumLikelihood <- function(data,LKmodel,n.parameters.LKmodel,...) {
  # recover()
  
  ## A boat-load of error checking
  
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
  
  # Make sure provided LKmodel is a function
  if ( !class(LKmodel) == "function" ) {
    stop("Argument LKmodel should be the name of a function (and not in quotes).
         For example, fitGeneralizedLotkaVolteraMaximumLikelihood(...,LKmodel=regularLK,...)")
  }
  
  # Make sure provided LKmodel is a function with n.parameters.LKmodel parameters in theta and 2 additional parameters
  if ( any(is.na(LKmodel(1,1,rep(1,n.parameters.LKmodel)))) ) {
    stop("Argument LKmodel should be a function that takes 2 + n.parameters.LKmodel arguments.\nThis example function has n.parameters.LKmodel = 2:
         \tmyExampleFunction <- function(prey.size,pred.size,theta){
         \t\tprey_slope <- exp(-theta[1]*prey.size)
         \t\tpred_slope <- exp(-theta[2]*pred.size)
         \t\treturn(c(prey_slope,pred_slope))
         \t}")
  }

  # Make sure KLmodel takes arguments with the correct names (this also catches errors that slip through the previous checks)
  if ( class(try(LKmodel(pred.size=1,prey.size=1,theta=rep(1,n.parameters.LKmodel)),silent=TRUE)) == "try-error" ) {
    stop("Argument LKmodel should be a function that takes 2 + n.parameters.LKmodel arguments.\nThe arguments must be named pred.size, prey.size, and theta. For example:
         \tmyExampleFunction <- function(prey.size,pred.size,theta){
         \t\tprey_slope <- exp(-theta[1]*prey.size)
         \t\tpred_slope <- exp(-theta[2]*pred.size)
         \t\treturn(c(prey_slope,pred_slope))
         \t}")
  }
  
  # Make sure provided LKmodel returns a numeric vector of length 2
  if ( length(LKmodel(pred.size=1,prey.size=1,theta=rep(1,n.parameters.LKmodel))) != 2 || class(LKmodel(pred.size=1,prey.size=1,theta=rep(1,n.parameters.LKmodel))) != "numeric" ) {
    stop("Argument LKmodel should be a function that returns a numeric vector of size 2")
  }
  
  ## First-pass estimation with derivatives, following Howard
  
  # The problem here is that starting with crazy optimization values can take forever and then yield crazy results
  # Part of this is due to the ability of the error model to soak up curve fitting failure
  # The solution is to start by using the observed derivatives to fit the non-error parameters
  # We can do this roughly and quickly with least squares
  # After this, we condition on these values and find reasonable starting values of the errors
  # Lastly, we take all the above as starting values and find the maximum likelihood solution
  
  # observed changes
  dH <- (data$prey[3:dim(data)[1]] - data$prey[1:(dim(data)[1]-2)])/(data$time[3:dim(data)[1]] - data$time[1:(dim(data)[1]-2)]) #nrow is slower than dim()[1] for esoteric reasons
  dP <- (data$pred[3:dim(data)[1]] - data$pred[1:(dim(data)[1]-2)])/(data$time[3:dim(data)[1]] - data$time[1:(dim(data)[1]-2)])
  
  observed_d_by_dt <- rbind(dH,dP)
  
  # recover()
  
  fitInitWithDerivs <- function(par) {
    predicted_d_by_dt <- sapply(2:(dim(data)[1]-1),function(i){
      LKmodel(prey.size=data$prey[i],pred.size=data$pred[i],theta=exp(par))
    })
    return(sum((predicted_d_by_dt - observed_d_by_dt)^2))
  }
  
  cat("Initializing model, please be patient!\n")
  # Now we initialize a handful of times, and find the best starting values
  # log(0.5) is a good starting guess for the simpler models, but can lead to bad optimization for the others
  init_attempt_n_tries_deterministic <- 10
  init_attempt_n_tries_random <- 50
  init_opt_attempts <- vector("list",init_attempt_n_tries_deterministic+init_attempt_n_tries_random)
  init_opt_scores <- numeric(init_attempt_n_tries_deterministic+init_attempt_n_tries_random)
  for (i in 1:init_attempt_n_tries_deterministic) {
    init_opt_attempts[[i]] <- optim(rep(log(i/10),n.parameters.LKmodel),fitInitWithDerivs)
    init_opt_scores[i] <- init_opt_attempts[[i]]$value
  }
  for (i in (1:init_attempt_n_tries_random)+init_attempt_n_tries_deterministic) {
    init_opt_attempts[[i]] <- optim(log(runif(n.parameters.LKmodel,0.1,1.1)),fitInitWithDerivs)
    init_opt_scores[i] <- init_opt_attempts[[i]]$value
  }
  
  # choose the best starting attempt
  best_init_index <- which(init_opt_scores == min(init_opt_scores))
  init_opts <- init_opt_attempts[[best_init_index]]
  
  ## Prepare for more thourough optimization
  
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
  
  # recover()
  
  ## Optimization step 2
  
  # Find area near peak for starting sizes and errors, conditioning on approximate parameters from least squares derivative fitting
  # If we don't do this, we can easily get lost in crazy parameters (optimization is unstable!)
  # A more principled approach is to start optimization many times in many places (parameter combinations), but this is too slow for a 2 hour lab, given the time to optimize
  fnLKInitAndVar <- function(par) {
    -logLikelihoodGeneralizedLK(data=data,
                                prey.initial.size=par[1],pred.initial.size=par[2],
                                sigma1=par[3],sigma2=par[4],
                                theta=init_opts$par,LKmodel=LKmodel,n.steps=n_steps_initial,parameters.are.log.scale=TRUE)
  }
  
  cat("Preparing to fit model, please be patient!\n")
  opts <- optim(c(log(data$prey[1]),log(data$pred[1]),log(0.3),log(0.3)),fnLKInitAndVar,method="BFGS")
  
  ## Optimization step 3
  
  fnLKforOptimML <- function(par) {
    -logLikelihoodGeneralizedLK(data=data,
                                prey.initial.size=par[1],pred.initial.size=par[2],
                                sigma1=par[3],sigma2=par[4],
                                theta=par[5:(4+n.parameters.LKmodel)],LKmodel=LKmodel,
                                n.steps=n_steps_final,parameters.are.log.scale=TRUE)
  }
  
  cat("Fitting model, please be patient!\n")
  opt_ml <- optim(c(opts$par,init_opts$par),fnLKforOptimML,method="BFGS",hessian=TRUE)
  
  ml_par <- exp(opt_ml$par) # move off log-scale
  names(ml_par) <- c("prey.initial.size","pred.initial.size","sigma1","sigma2",paste0("theta",1:n.parameters.LKmodel))
  
  # recover()
  
  # Compute confidence intervals (keeping in mind we worked on log-scale for all optimization)
  asymptotic_sd <- sqrt(1/diag(opt_ml$hessian))
  asymptotic_025 <- qlnorm(0.025,opt_ml$par,asymptotic_sd)
  asymptotic_975 <- qlnorm(0.975,opt_ml$par,asymptotic_sd)

  asymptotic_ci <- cbind(asymptotic_025,asymptotic_975)
  colnames(asymptotic_ci) <- c("2.5%","97.5%")
  rownames(asymptotic_ci) <- c("prey.initial.size","pred.initial.size","sigma1","sigma2",paste0("theta",1:n.parameters.LKmodel))
  
  return(list(logLikelihood=-opt_ml$value,maximum.likelihood.parameter.estimates=ml_par,confidence.intervals=asymptotic_ci))
  
}

compareModels <- function(fitted.models.list) {
  # recover()
  
  ## Check inputs
  # Input of correct class?
  if ( class(fitted.models.list) != "list" ) {
    stop("compareModels requires input is a list")
  }
  # Input has likelihood values?
  all_input_names <- attributes(unlist(fitted.models.list))$names
  if ( length(grep("logLikelihood",all_input_names)) != length(fitted.models.list) ) {
    stop("compareModels requires the full output of fitGeneralizedLotkaVolteraMaximumLikelihood for each model")
  }
  # Input has believable structure? Should be twice as many confidence intervals as estimated parameters
  if ( length(grep("confidence.intervals",all_input_names)) != 2*length(grep("maximum.likelihood.parameter.estimates",all_input_names))) {
    stop("compareModels requires the full output of fitGeneralizedLotkaVolteraMaximumLikelihood for each model")
  }
  
  n_models <- length(fitted.models.list)
  
  ## Calculate AIC scores
  aic_scores <- numeric(n_models)
  for (i in 1:n_models) {
    aic_scores[i] <- -2*fitted.models.list[[i]]$logLikelihood + 2*length(fitted.models.list[[i]]$maximum.likelihood.parameter.estimates)
  }
  ranks <- order(aic_scores) # keep order so we can name results later
  aic_scores <- sort(aic_scores) # sort for easy comparison
  
  ## AIC relative likelihoods
  aic_rl <- exp((aic_scores[1] - aic_scores)/2) # first is the best/minimum
  names(aic_rl) <- paste0("Model ",ranks)
  
  cat("Relative likelihoods of all models under consideration (compared to best):\n")
  print(aic_rl)
}



assessGeneralizedLKModelFit <- function(data,LKmodel,fitted.model,n.sim=100,...) {
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
  n_time_points <- length(data$time)
  
  # Pick a time step size
  # This default will work reasonably well for data with two peaks
  # We want about 4000-5000 steps
  total_time <- max(data$time) - min(data$time)
  
  # Step size for trajectory simulation
  steps_per_interval <- floor(5000/total_time)
  n_steps <- total_time*steps_per_interval
  
  # Get parameters in a form we can use
  is_theta <- grepl("theta",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_prey_init <- grepl("prey.initial.size",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_pred_init <- grepl("pred.initial.size",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_sigma1 <- grepl("sigma1",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_sigma2 <- grepl("sigma2",names(fitted.model$maximum.likelihood.parameter.estimates))
  fitted_theta <- fitted.model$maximum.likelihood.parameter.estimates[is_theta]
  fitted_prey_init <- fitted.model$maximum.likelihood.parameter.estimates[is_prey_init]
  fitted_pred_init <- fitted.model$maximum.likelihood.parameter.estimates[is_pred_init]
  fitted_sigma1 <- fitted.model$maximum.likelihood.parameter.estimates[is_sigma1]
  fitted_sigma2 <- fitted.model$maximum.likelihood.parameter.estimates[is_sigma2]
  
  # Simulate trajectories from model
  pred_sim <- matrix(NA,nrow=n_time_points,ncol=n.sim)
  prey_sim <- matrix(NA,nrow=n_time_points,ncol=n.sim)
  pb <- txtProgressBar(min = 0, max = n.sim, style = 3) # report progress
  for (i in 1:n.sim) {
    # Simulate true trajectory
    sim_traj <- sampleGeneralizedLKTrajectory(prey.initial.size=fitted_prey_init,pred.initial.size=fitted_pred_init,theta=fitted_theta,LKmodel=LKmodel,sampled.times=data$time,n.steps=n_steps)

    # add measurement error
    prey_sim[,i] <- rlnorm(n_time_points,log(sim_traj$prey_trajectory),fitted_sigma1)
    pred_sim[,i] <- rlnorm(n_time_points,log(sim_traj$pred_trajectory),fitted_sigma2)
    setTxtProgressBar(pb, i) # report progress
  }
  
  # Plot model

  # # so we can put axes on same scale
  # plot_max <- max(c(max(prey_sim),max(pred_sim))) # since there's measurement error, it's possible that there's a super-high predator pop size
  
  par(mfrow=c(2,1),mai=c(0.8,0.8,0.2,0.2),xpd=TRUE)
  matplot(x=data$time,prey_sim,col=prey_color,type="l",lty=1,xlab="time",ylab="prey pop size")
  lines(data$time,data$prey,lwd=3)

  matplot(x=data$time,pred_sim,col=pred_color,type="l",lty=1,xlab="time",ylab="predator pop size")
  lines(data$time,data$pred,lwd=3)
  par(mfrow=c(1,1))
  
}
