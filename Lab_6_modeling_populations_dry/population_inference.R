# Takes the following arguments:
#      data: the data, as a dataframe with columns named time(s), pred(ator), and prey
plotPredatorPrey <- function(data) {
  par(mfrow=c(2,1),mai=c(0.8,0.8,0.2,0.2),xpd=TRUE)
  plot(data$time,data$prey,type="l",lwd=2,xlab="time",ylab="prey pop size")
  
  plot(data$time,data$pred,type="l",lwd=2,xlab="time",ylab="predator pop size")
  par(mfrow=c(1,1))
}

# This function fits population models with maximum likelihood
# Models are specified with 4 functions, 1 each to describe birth and death rates of predators and prey
# These functions take as input the predator and prey population size and a vector of at least one other parameter
# Takes the following arguments:
#      data: the data, as a dataframe with columns named time(s), pred(ator), and prey
#      preyBirthFxn: a function that returns the prey birth rate for given prey and predator population sizes and at least one additional parameter
#      preyDeathFxn: a function that returns the prey death rate for given prey and predator population sizes and at least one additional parameter
#      predBirthFxn: a function that returns the predator birth rate for given prey and predator population sizes and at least one additional parameter
#      predDeathFxn: a function that returns the predator death rate for given prey and predator population sizes and at least one additional parameter
#      n.parameters.in.functions: a vector of length 4, the number of parameters for the functions in order preyBirthFxn,preyDeathFxn,predBirthFxn,predDeathFxn
#      n.initialization.attempts: how many starting points to consider (default 2000, some datasets are unstable and need more)
# Optionally, takes
#      initialize.grid.upper:
#      initialize.grid.lower:
fitGeneralizedLotkaVolteraMaximumLikelihood <- function(data,
                                                        preyBirthFxn,
                                                        preyDeathFxn,
                                                        predBirthFxn,
                                                        predDeathFxn,
                                                        n.parameters.in.functions,
                                                        n.initialization.attempts=2000,
                                                        ...) {
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
  
  # There are many things we have to check for the input functions, so we outsource it to a different function
  # The error-checking function will fail loudly and stop everything from proceeding
  assertAllFunctionsPassAllTests(preyBirthFxn,preyDeathFxn,predBirthFxn,predDeathFxn,n.parameters.in.functions)
  
  ## First-pass estimation with derivatives, following Howard
  
  # The problem here is that starting with crazy optimization values can take forever and then yield crazy results
  # Part of this is due to the ability of the error model to soak up curve fitting failure
  # The solution is to start by using the observed derivatives to fit the non-error parameters
  # We can do this roughly and quickly with least squares
  # After this, we condition on these values and find reasonable starting values of the errors
  # Lastly, we take all the above as starting values and find the maximum likelihood solution
  
  # The initialization checks a grid of parameter values, spaced in powers of 10, from 10^initialize_grid_lower to 10^initialize_grid_upper
  # The default is a lower of -5 and upper of 5
  # For most of these models that's fine, for some it may help to expand the upper bound (especially for large population sizes when estimating carrying capacities)
  if ( hasArg(initialize.grid.upper) ) {
    initialize_grid_upper <- list(...)$initialize.grid.upper
  } else {
    initialize_grid_upper <- 5
  }
  if (hasArg(initialize.grid.lower)) {
    initialize_grid_lower <- list(...)$initialize.grid.lower
  } else {
    initialize_grid_lower <- -5
  }
  cat("Initializing model, please be patient!\n")
  
  init_opts <- initializeLVGeneral(data,
                                   preyBirthFxn,
                                   preyDeathFxn,
                                   predBirthFxn,
                                   predDeathFxn,
                                   n.parameters.in.functions,
                                   n.initialization.attempts,
                                   10,
                                   initialize_grid_upper,
                                   initialize_grid_lower)
  
  # recover()
  
  ## Prepare for more thourough optimization
  
  # Pick a time step size
  # This default will work reasonably well for data with two peaks
  # We want about 4000-5000 steps
  total_time <- max(data$time) - min(data$time)
  
  # For final optimization, smaller steps
  steps_per_interval_final <- floor(4000/total_time)
  n_steps_final <- total_time*steps_per_interval_final
  
  # For initial, less precise is fine
  steps_per_interval_initial <- floor(1000/total_time)
  n_steps_initial <- total_time*steps_per_interval_initial
  
  # recover()
  
  ## Optimization step 2
  
  # Find area near peak for starting sizes and errors, conditioning on approximate parameters from least squares derivative fitting
  # If we don't do this, we can easily get lost in crazy parameters (optimization is unstable!)
  # A more principled approach is to start optimization many times in many places (parameter combinations), but this is too slow for a 2 hour lab, given the time to optimize
  cat("Preparing to fit model, please be patient!\n")
  init_opts_2 <- initializeLVPopsAndErrors(data,
                                           preyBirthFxn,
                                           preyDeathFxn,
                                           predBirthFxn,
                                           predDeathFxn,
                                           n.parameters.in.functions,
                                           n_steps_initial,
                                           init_opts)
  
  ## Optimization step 3
  
  # Establish where the parameters will be during optimization
  prey_birth_indices <- 4+1:(0+n.parameters.in.functions[1])
  prey_death_indices <- 4+(1+n.parameters.in.functions[1]):(0+sum(n.parameters.in.functions[1:2]))
  pred_birth_indices <- 4+(1+sum(n.parameters.in.functions[1:2])):(0+sum(n.parameters.in.functions[1:3]))
  pred_death_indices <- 4+(1+sum(n.parameters.in.functions[1:3])):(0+sum(n.parameters.in.functions[1:4]))
  
  fnLKforOptimML <- function(par) {
    -logLikelihoodGeneralizedLV(data=data,
                                prey.initial.size=par[1],
                                pred.initial.size=par[2],
                                sigma1=par[3],
                                sigma2=par[4],
                                preyBirthFxn=preyBirthFxn,
                                preyDeathFxn=preyDeathFxn,
                                predBirthFxn=predBirthFxn,
                                predDeathFxn=predDeathFxn,
                                prey.birth.parameters=par[prey_birth_indices],
                                prey.death.parameters=par[prey_death_indices],
                                pred.birth.parameters=par[pred_birth_indices],
                                pred.death.parameters=par[pred_death_indices],
                                n.steps=n_steps_final,
                                parameters.are.log.scale=TRUE)
  }
  
  cat("Fitting model, please be patient!\n")
  opt_ml <- optim(c(init_opts_2,init_opts),fnLKforOptimML,method="BFGS",hessian=TRUE)
  
  # recover()
  
  ml_par <- exp(opt_ml$par) # move off log-scale
  names(ml_par) <- c("prey.initial.size",
                     "pred.initial.size",
                     "sigma1",
                     "sigma2",
                     paste0("prey.birth.parameters",1:n.parameters.in.functions[1]),
                     paste0("prey.death.parameters",1:n.parameters.in.functions[2]),
                     paste0("pred.birth.parameters",1:n.parameters.in.functions[3]),
                     paste0("pred.death.parameters",1:n.parameters.in.functions[4]))
  
  # Compute confidence intervals (keeping in mind we worked on log-scale for all optimization)
  asymptotic_sd <- sqrt(1/diag(opt_ml$hessian))
  asymptotic_025 <- qlnorm(0.025,opt_ml$par,asymptotic_sd)
  asymptotic_975 <- qlnorm(0.975,opt_ml$par,asymptotic_sd)
  
  asymptotic_ci <- cbind(asymptotic_025,asymptotic_975)
  colnames(asymptotic_ci) <- c("2.5%","97.5%")
  rownames(asymptotic_ci) <- names(ml_par)
  
  return(list(logLikelihood=-opt_ml$value,maximum.likelihood.parameter.estimates=ml_par,confidence.intervals=asymptotic_ci))
  
}

# This is a helper function for the full ML analysis
# It uses the observed derivatives and least squares to fit the birth and death rate parameters
# These values are not the ones ML will return, but they are computable much more quickly
# Since this setup takes the population sizes at a given time to be true (no error) it also allows us to estimate predator and prey separately, a further speed boost
# Thus, we can explore a wider space of staring values here, and then polish with ML (which is much slower here)
# Takes the following arguments:
#      data: the data, as a dataframe with columns named time(s), pred(ator), and prey
#      preyBirthFxn: a function that returns the prey birth rate for given prey and predator population sizes and at least one additional parameter
#      preyDeathFxn: a function that returns the prey death rate for given prey and predator population sizes and at least one additional parameter
#      predBirthFxn: a function that returns the predator birth rate for given prey and predator population sizes and at least one additional parameter
#      predDeathFxn: a function that returns the predator death rate for given prey and predator population sizes and at least one additional parameter
#      n.parameters.in.functions: a vector of length 4, the number of parameters for the functions in order preyBirthFxn,preyDeathFxn,predBirthFxn,predDeathFxn
#      n.starting.points: how many random points should we start considering (we will take only a few for further consideration)?
#      n.starting.points.to.optimize: How many of the above should we consider in more depth? We will take the n best
# Returns a vector of starting parameters in the order prey birth, prey death, pred birth, pred death 
initializeLVGeneral <- function(data,
                                preyBirthFxn,
                                preyDeathFxn,
                                predBirthFxn,
                                predDeathFxn,
                                n.parameters.in.functions,
                                n.starting.points,
                                n.starting.points.to.optimize,
                                initialize.grid.upper=5,
                                initialize.grid.lower=-5) {
  
  # recover()
  
  # observed changes at given times
  observed_dPrey <- (data$prey[3:dim(data)[1]] - data$prey[1:(dim(data)[1]-2)])/(data$time[3:dim(data)[1]] - data$time[1:(dim(data)[1]-2)]) #nrow is slower than dim()[1] for esoteric reasons
  observed_dPred <- (data$pred[3:dim(data)[1]] - data$pred[1:(dim(data)[1]-2)])/(data$time[3:dim(data)[1]] - data$time[1:(dim(data)[1]-2)])
  
  ## Least squares for prey trajectory
  
  # Establish where the parameters will be during optimization for prey parameters
  prey_birth_indices <- 1:(0+n.parameters.in.functions[1])
  prey_death_indices <- (1+n.parameters.in.functions[1]):(0+sum(n.parameters.in.functions[1:2]))
  
  # recover()
  
  # Calculates the sum of squared errors from model/predicted derivatives to observed derivatives
  sumOfSquaresDerivativesPrey <- function(par) {
    # calculate predicted slopes given model parameters
    predicted_dPrey <- sapply(2:(dim(data)[1]-1),function(i){
      getPreySlope(prey.size=data$prey[i],
                   pred.size=data$pred[i],
                   preyBirthFxn=preyBirthFxn,
                   preyDeathFxn=preyDeathFxn,
                   prey.birth.parameters=par[prey_birth_indices],
                   prey.death.parameters=par[prey_death_indices])
    })

    # Return sum of squared errors (from model predictions at these parameters to observed values)
    return(sum((predicted_dPrey - observed_dPrey)^2))
  }
  
  # Layout points for initializing (we may not hit the exact number the user asked for)
  # Space points logarithmically in a reasonable range
  n_prey_param <- sum(n.parameters.in.functions[1:2])
  n_starting_per_prey_param <- (n.starting.points)^(1/n_prey_param)
  vals <- exp(seq(initialize.grid.lower,initialize.grid.upper,length.out=n_starting_per_prey_param))
  prey_start_list <- vector("list",n_prey_param)
  for (i in 1:n_prey_param) {
    prey_start_list[[i]] <- vals
  }
  prey_grid <- expand.grid(prey_start_list)
  
  # Loop over rows, calculate how decent a fit this is
  sum_squares_prey <- apply(prey_grid,1,sumOfSquaresDerivativesPrey)
  
  prey_starting_points_rankings <- order(sum_squares_prey)
  
  # New function for optimizing (have to put params on log scale so they aren't bounded)
  sumOfSquaresDerivativesPrey <- function(par) {
    
    par <- exp(par)
    
    # calculate predicted slopes given model parameters
    predicted_dPrey <- sapply(2:(dim(data)[1]-1),function(i){
      getPreySlope(prey.size=data$prey[i],
                   pred.size=data$pred[i],
                   preyBirthFxn=preyBirthFxn,
                   preyDeathFxn=preyDeathFxn,
                   prey.birth.parameters=par[prey_birth_indices],
                   prey.death.parameters=par[prey_death_indices])
    })
    
    # Return sum of squared errors (from model predictions at these parameters to observed values)
    return(sum((predicted_dPrey - observed_dPrey)^2))
  }
  
  # loop over the best set of potential starting points, and optimize
  prey_init_opts <- lapply(prey_starting_points_rankings[1:n.starting.points.to.optimize],function(i){
    try(optim(log(prey_grid[i,]),sumOfSquaresDerivativesPrey,control=list(reltol=1e-4)),silent=TRUE)
  })
  
  # Find best of these, re-optimize with smaller tolerance
  # Then these are our initial prey optimization parameters
  sum_squares_prey <- numeric(n.starting.points.to.optimize)
  for (i in 1:n.starting.points.to.optimize) {
    sum_squares_prey[i] <- ifelse(class(prey_init_opts[[i]]) == "try-error", Inf, prey_init_opts[[i]]$value)
  }
  prey_start <- prey_init_opts[[which(sum_squares_prey == min(sum_squares_prey))[1]]]$par # We take [1] of the which() statement in case there are multiple values with the same endpoint (they're equally good, take any!)
  
  prey_init_opts <- optim(prey_start,sumOfSquaresDerivativesPrey)
  prey_start <- prey_init_opts$par
  
  ## Least squares for predator trajectory
  
  # Establish where the parameters will be during optimization for predator parameters
  pred_birth_indices <- 1:(0+n.parameters.in.functions[3])
  pred_death_indices <- (1+n.parameters.in.functions[3]):(0+sum(n.parameters.in.functions[3:4]))
  
  sumOfSquaresDerivativesPred <- function(par) {
    # calculate predicted slopes given model parameters
    predicted_dPred <- sapply(2:(dim(data)[1]-1),function(i){
      getPredSlope(prey.size=data$prey[i],
                   pred.size=data$pred[i],
                   predBirthFxn=predBirthFxn,
                   predDeathFxn=predDeathFxn,
                   pred.birth.parameters=par[pred_birth_indices],
                   pred.death.parameters=par[pred_death_indices])
    })
    
    # Return sum of squared errors (from model predictions at these parameters to observed values)
    return(sum((predicted_dPred - observed_dPred)^2))
  }
  
  # Layout points for initializing (we may not hit the exact number the user asked for)
  # Space points logarithmically in a reasonable range
  n_pred_param <- sum(n.parameters.in.functions[3:4])
  n_starting_per_pred_param <- (n.starting.points)^(1/n_pred_param)
  vals <- exp(seq(initialize.grid.lower,initialize.grid.upper,length.out=n_starting_per_pred_param))
  pred_start_list <- vector("list",n_pred_param)
  for (i in 1:n_pred_param) {
    pred_start_list[[i]] <- vals
  }
  pred_grid <- expand.grid(pred_start_list)
  
  # Loop over rows, calculate how decent a fit this is
  sum_squares_pred <- apply(pred_grid,1,sumOfSquaresDerivativesPred)
  
  pred_starting_points_rankings <- order(sum_squares_pred)
  
  # New function for optimizing (have to put params on log scale so they aren't bounded)
  
  sumOfSquaresDerivativesPred <- function(par) {
    # Get off log scale
    par <- exp(par)
    
    # calculate predicted slopes given model parameters
    predicted_dPred <- sapply(2:(dim(data)[1]-1),function(i){
      getPredSlope(prey.size=data$prey[i],
                   pred.size=data$pred[i],
                   predBirthFxn=predBirthFxn,
                   predDeathFxn=predDeathFxn,
                   pred.birth.parameters=par[pred_birth_indices],
                   pred.death.parameters=par[pred_death_indices])
    })
    
    # Return sum of squared errors (from model predictions at these parameters to observed values)
    return(sum((predicted_dPred - observed_dPred)^2))
  }
  
  # loop over the best set of potential starting points, and optimize
  pred_init_opts <- lapply(pred_starting_points_rankings[1:n.starting.points.to.optimize],function(i){
    try(optim(log(pred_grid[i,]),sumOfSquaresDerivativesPred,control=list(reltol=1e-4)),silent=TRUE)
  })
  
  # Find best of these, those are our initial pred optimization parameters
  sum_squares_pred <- numeric(n.starting.points.to.optimize)
  for (i in 1:n.starting.points.to.optimize) {
    sum_squares_pred[i] <- ifelse(class(pred_init_opts[[i]]) == "try-error", Inf, pred_init_opts[[i]]$value)
  }
  pred_start <- pred_init_opts[[which(sum_squares_pred == min(sum_squares_pred))[1]]]$par# We take [1] of the which() statement in case there are multiple values with the same endpoint (they're equally good, take any!)
  
  pred_init_opts <- optim(log(pred_grid[i,]),sumOfSquaresDerivativesPred)
  pred_start <- pred_init_opts$par
  
  return(c(prey_start,pred_start))
  
}

# This is a helper function for the full ML analysis
# It initializes reasonable values of error terms and initial population sizes before the full ML analysis
# Takes the following arguments:
#      data: the data, as a dataframe with columns named time(s), pred(ator), and prey
#      preyBirthFxn: a function that returns the prey birth rate for given prey and predator population sizes and at least one additional parameter
#      preyDeathFxn: a function that returns the prey death rate for given prey and predator population sizes and at least one additional parameter
#      predBirthFxn: a function that returns the predator birth rate for given prey and predator population sizes and at least one additional parameter
#      predDeathFxn: a function that returns the predator death rate for given prey and predator population sizes and at least one additional parameter
#      n.parameters.in.functions: a vector of length 4, the number of parameters for the functions in order preyBirthFxn,preyDeathFxn,predBirthFxn,predDeathFxn
#      n.starting.points: how many random points should we start considering (we will take only a few for further consideration)?
#      n.starting.points.to.optimize: How many of the above should we consider in more depth? We will take the n best
# Returns a vector of starting parameters in the order prey birth, prey death, pred birth, pred death 
initializeLVPopsAndErrors <- function(data,
                                      preyBirthFxn,
                                      preyDeathFxn,
                                      predBirthFxn,
                                      predDeathFxn,
                                      n.parameters.in.functions,
                                      n_steps_initial,
                                      init_opts) {

  # Establish where the parameters will be during optimization
  prey_birth_indices <- 1:(0+n.parameters.in.functions[1])
  prey_death_indices <- (1+n.parameters.in.functions[1]):(0+sum(n.parameters.in.functions[1:2]))
  pred_birth_indices <- (1+sum(n.parameters.in.functions[1:2])):(0+sum(n.parameters.in.functions[1:3]))
  pred_death_indices <- (1+sum(n.parameters.in.functions[1:3])):(0+sum(n.parameters.in.functions[1:4]))
  
  
  fnInitPops <- function(par) {
    -logLikelihoodGeneralizedLV(data=data,
                                prey.initial.size=par[1],
                                pred.initial.size=par[2],
                                sigma1=-2.3, #~= log(0.1)
                                sigma2=-2.3, #~= log(0.1)
                                preyBirthFxn=preyBirthFxn,
                                preyDeathFxn=preyDeathFxn,
                                predBirthFxn=predBirthFxn,
                                predDeathFxn=predDeathFxn,
                                prey.birth.parameters=init_opts[prey_birth_indices],
                                prey.death.parameters=init_opts[prey_death_indices],
                                pred.birth.parameters=init_opts[pred_birth_indices],
                                pred.death.parameters=init_opts[pred_death_indices],
                                n.steps=n_steps_initial,
                                parameters.are.log.scale=TRUE)
  }
  
  opt_pops <- optim(c(log(data$prey[1]),log(data$pred[1])),fnInitPops,method="BFGS")
  
  # Grab the path laid out by these parameters, will use to get ML error terms
  path <- sampleGeneralizedLKTrajectory(prey.initial.size=exp(opt_pops$par[1]),
                                        pred.initial.size=exp(opt_pops$par[2]),
                                        preyBirthFxn=preyBirthFxn,
                                        preyDeathFxn=preyDeathFxn,
                                        predBirthFxn=predBirthFxn,
                                        predDeathFxn=predDeathFxn,
                                        prey.birth.parameters=exp(init_opts[prey_birth_indices]),
                                        prey.death.parameters=exp(init_opts[prey_death_indices]),
                                        pred.birth.parameters=exp(init_opts[pred_birth_indices]),
                                        pred.death.parameters=exp(init_opts[pred_death_indices]),
                                        sampled.times=data$time,
                                        n.steps=n_steps_initial)
  
  # Optimize first error term (for prey path)
  fnInitErrors1 <- function(par) {
    -sum(dlnorm(data$prey,log(path$prey_trajectory),par,log=TRUE))
  }
  
  opt_error1 <- optimize(fnInitErrors1,interval=c(0,3))
  
  # Optimize second error term (for predator path)
  fnInitErrors2 <- function(par) {
    -sum(dlnorm(data$pred,log(path$pred_trajectory),par,log=TRUE))
  }
  
  opt_error2 <- optimize(fnInitErrors2,interval=c(0,3))
  
  # Clean up by jointly optimizing all 4 things, starting at current best
  fnInitPopsAndErrors <- function(par) {
    # cat(exp(par),"\n")
    -logLikelihoodGeneralizedLV(data=data,
                                prey.initial.size=par[1],
                                pred.initial.size=par[2],
                                sigma1=par[3],
                                sigma2=par[4],
                                preyBirthFxn=preyBirthFxn,
                                preyDeathFxn=preyDeathFxn,
                                predBirthFxn=predBirthFxn,
                                predDeathFxn=predDeathFxn,
                                prey.birth.parameters=init_opts[prey_birth_indices],
                                prey.death.parameters=init_opts[prey_death_indices],
                                pred.birth.parameters=init_opts[pred_birth_indices],
                                pred.death.parameters=init_opts[pred_death_indices],
                                n.steps=n_steps_initial,
                                parameters.are.log.scale=TRUE)
  }
  
  opt_both <- optim(c(opt_pops$par,log(opt_error1$minimum),log(opt_error2$minimum)),fnInitPopsAndErrors,method="BFGS")
  
  return(opt_both$par)
  
}
# Calculates the likelihood for a generalized class of Lotka-Voltera models
# Takes the following argumens:
#      data: the data, as a dataframe with columns named time(s), pred(ator), and prey
#      prey.initial.size: the population size for the prey at the first time point
#      pred.initial.size: the population size for the predator at the first time point
#      preyBirthFxn: a function that returns the prey birth rate for given prey and predator population sizes and at least one additional parameter
#      preyDeathFxn: a function that returns the prey death rate for given prey and predator population sizes and at least one additional parameter
#      predBirthFxn: a function that returns the predator birth rate for given prey and predator population sizes and at least one additional parameter
#      predDeathFxn: a function that returns the predator death rate for given prey and predator population sizes and at least one additional parameter
#      prey.birth.parameters: the parameters for the model of prey birth rate
#      prey.death.parameters: the parameters for the model of prey death rate
#      pred.birth.parameters: the parameters for the model of predator birth rate
#      pred.death.parameters: the parameters for the model of predator death rate
#      sigma1: the value of sigma for the (lognormal) error of prey measurements
#      sigma2: the value of sigma for the (lognormal) error of predator measurements
#      parameters.are.log.scale: are the parameter values given on the log scale? (They will be for optimization)
logLikelihoodGeneralizedLV <- function(data,
                                       prey.initial.size,
                                       pred.initial.size,
                                       preyBirthFxn,
                                       preyDeathFxn,
                                       predBirthFxn,
                                       predDeathFxn,
                                       prey.birth.parameters,
                                       prey.death.parameters,
                                       pred.birth.parameters,
                                       pred.death.parameters,
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
    prey.birth.parameters <- exp(prey.birth.parameters)
    prey.death.parameters <- exp(prey.death.parameters)
    pred.birth.parameters <- exp(pred.birth.parameters)
    pred.death.parameters <- exp(pred.death.parameters)
    sigma1 <- exp(sigma1)
    sigma2 <- exp(sigma2)
  }
  
  # Need to know what we think the pop sizes were at the times we sampled
  path <- sampleGeneralizedLKTrajectory(prey.initial.size=prey.initial.size,
                                        pred.initial.size=pred.initial.size,
                                        preyBirthFxn=preyBirthFxn,
                                        preyDeathFxn=preyDeathFxn,
                                        predBirthFxn=predBirthFxn,
                                        predDeathFxn=predDeathFxn,
                                        prey.birth.parameters=prey.birth.parameters,
                                        prey.death.parameters=prey.death.parameters,
                                        pred.birth.parameters=pred.birth.parameters,
                                        pred.death.parameters=pred.death.parameters,
                                        sampled.times=data$time,
                                        n.steps=n.steps)
  
  # Sometimes the path will have a population crash to 0, but it keeps slogging along and creates NaNs
  # Here, if a population goes to 0, the probability of observing >0 individuals is 0 (lnL = -Inf)
  # So we catch this and return -<very big number> (as a hack to avoid non-finite differences in optimization)
  flat_path <- unlist(path) # path as a flat vector, can't use is.nan() and related functions on data.frames
  if ( sum(is.nan(flat_path)) > 0 || sum(is.infinite(flat_path)) > 0 || sum(flat_path == 0) > 0 ) {
    return(-.Machine$integer.max)
  }
  
  # Now it's just a measurement model
  lnL <- sum(dlnorm(data$prey,log(path$prey_trajectory),sigma1,log=TRUE)) + # measurements of prey pop sizes
    sum(dlnorm(data$pred,log(path$pred_trajectory),sigma2,log=TRUE)) # measurements of pred pop sizes
  
  return(lnL)
  
}

# Here we simulate a trajectory for the predator and prey population sizes using models for birth and death rates of predator and prey populations
# In other words, this function can plot the trajectory for just about any (deterministic) model of predator and prey population cycling that depends only on predator and prey population sizes
# Takes the following argumens:
#      prey.initial.size: the population size for the prey at the first time point
#      pred.initial.size: the population size for the predator at the first time point
#      preyBirthFxn: a function that returns the prey birth rate for given prey and predator population sizes and at least one additional parameter
#      preyDeathFxn: a function that returns the prey death rate for given prey and predator population sizes and at least one additional parameter
#      predBirthFxn: a function that returns the predator birth rate for given prey and predator population sizes and at least one additional parameter
#      predDeathFxn: a function that returns the predator death rate for given prey and predator population sizes and at least one additional parameter
#      prey.birth.parameters: the parameters for the model of prey birth rate
#      prey.death.parameters: the parameters for the model of prey death rate
#      pred.birth.parameters: the parameters for the model of predator birth rate
#      pred.death.parameters: the parameters for the model of predator death rate
sampleGeneralizedLKTrajectory <- function(prey.initial.size,
                                          pred.initial.size,
                                          preyBirthFxn,
                                          preyDeathFxn,
                                          predBirthFxn,
                                          predDeathFxn,
                                          prey.birth.parameters,
                                          prey.death.parameters,
                                          pred.birth.parameters,
                                          pred.death.parameters,
                                          sampled.times,
                                          n.steps) {
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
    prey_slope <- getPreySlope(prey.size=prey_trajectory[t],
                               pred.size=pred_trajectory[t],
                               preyBirthFxn=preyBirthFxn,
                               preyDeathFxn=preyDeathFxn,
                               prey.birth.parameters=prey.birth.parameters,
                               prey.death.parameters=prey.death.parameters)
    pred_slope <- getPredSlope(prey.size=prey_trajectory[t],
                               pred.size=pred_trajectory[t],
                               predBirthFxn=predBirthFxn,
                               predDeathFxn=predDeathFxn,
                               pred.birth.parameters=pred.birth.parameters,
                               pred.death.parameters=pred.death.parameters)
    prey_trajectory[t+1] <- prey_trajectory[t]+prey_slope*h
    pred_trajectory[t+1] <- pred_trajectory[t]+pred_slope*h
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

# This is a helper function that combines birth and death rates into an overall growth rate
getPreySlope <- function(prey.size,pred.size,preyBirthFxn,preyDeathFxn,prey.birth.parameters,prey.death.parameters) {
  per_capita_slope <- preyBirthFxn(pred.size=pred.size,prey.size=prey.size,prey.birth.parameters=prey.birth.parameters) - 
    preyDeathFxn(pred.size=pred.size,prey.size=prey.size,prey.death.parameters=prey.death.parameters)
  return(per_capita_slope * prey.size)
}

# This is a helper function that combines birth and death rates into an overall growth rate
getPredSlope <- function(prey.size,pred.size,predBirthFxn,predDeathFxn,pred.birth.parameters,pred.death.parameters) {
  per_capita_slope <- predBirthFxn(pred.size=pred.size,prey.size=prey.size,pred.birth.parameters=pred.birth.parameters) - 
    predDeathFxn(pred.size=pred.size,prey.size=prey.size,pred.death.parameters=pred.death.parameters)
  return(per_capita_slope * pred.size)
}

# This function takes in fitted models and uses AIC (the Akaike Information Criterion) to compare models
# The AIC seeks to balance the simplicity of a model (fewer parameters) with explanatory power (high likelihood)
# This function returns the relative likelihoods of the input models, it does not have a particularly useful interpretation, other than as a plausibility relative to the best fit model
# Technically speaking, the relative likelihood can be interpreted as being proportional to the probability that the ith model minimizes the (estimated) information loss
# Takes the following argumens:
#      fitted.models.list: a list of the fitted models from fitGeneralizedLotkaVolteraMaximumLikelihood (can make as list(fitted1,fitted2))
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

# This function 
# Takes the following argumens:
#      prey.initial.size: the population size for the prey at the first time point
#      pred.initial.size: the population size for the predator at the first time point
#      preyBirthFxn: a function that returns the prey birth rate for given prey and predator population sizes and at least one additional parameter
#      preyDeathFxn: a function that returns the prey death rate for given prey and predator population sizes and at least one additional parameter
#      predBirthFxn: a function that returns the predator birth rate for given prey and predator population sizes and at least one additional parameter
#      predDeathFxn: a function that returns the predator death rate for given prey and predator population sizes and at least one additional parameter
#      fitted.models.list: a list of the fitted models from fitGeneralizedLotkaVolteraMaximumLikelihood (can make as list(fitted1,fitted2))
#      n.sim: the number of trajectories to simulate
# Optionally, takes
#      pred.color: a hexadecimal color used to plot simulated trajectories
#      prey.color: a hexadecimal color used to plot simulated trajectories
#      relative.y.limit: 
visualizeGeneralizedLVModelFit <- function(data,
                                           preyBirthFxn,
                                           preyDeathFxn,
                                           predBirthFxn,
                                           predDeathFxn,
                                           fitted.model,
                                           n.sim=100,
                                           ...) {
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
  
  ## Check inputs
  if ( class(data) != "data.frame" ) {
    stop("Error in visualizeGeneralizedLVModelFit, data must be a data.frame")
  }
  if  ( class(fitted.model) != "list" ) {
    stop("Error in visualizeGeneralizedLVModelFit, is fitted.model the output of fitGeneralizedLotkaVolteraMaximumLikelihood?")
  }
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
  is_prey_birth <- grepl("prey.birth.parameters",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_prey_death <- grepl("prey.death.parameters",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_pred_birth <- grepl("pred.birth.parameters",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_pred_death <- grepl("pred.death.parameters",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_prey_init <- grepl("prey.initial.size",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_pred_init <- grepl("pred.initial.size",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_sigma1 <- grepl("sigma1",names(fitted.model$maximum.likelihood.parameter.estimates))
  is_sigma2 <- grepl("sigma2",names(fitted.model$maximum.likelihood.parameter.estimates))
  fitted_prey_birth <- fitted.model$maximum.likelihood.parameter.estimates[is_prey_birth]
  fitted_prey_death <- fitted.model$maximum.likelihood.parameter.estimates[is_prey_death]
  fitted_pred_birth <- fitted.model$maximum.likelihood.parameter.estimates[is_pred_birth]
  fitted_pred_death <- fitted.model$maximum.likelihood.parameter.estimates[is_pred_death]
  fitted_prey_init <- fitted.model$maximum.likelihood.parameter.estimates[is_prey_init]
  fitted_pred_init <- fitted.model$maximum.likelihood.parameter.estimates[is_pred_init]
  fitted_sigma1 <- fitted.model$maximum.likelihood.parameter.estimates[is_sigma1]
  fitted_sigma2 <- fitted.model$maximum.likelihood.parameter.estimates[is_sigma2]
  
  # Check we can find birth-death parameters
  if ( class(fitted_prey_birth) != "numeric" || class(fitted_prey_death) != "numeric" || class(fitted_pred_birth) != "numeric" || class(fitted_pred_death) != "numeric") {
    stop("Error in visualizeGeneralizedLVModelFit, cannot find birth-death parameters in fitted.model, is this argument the ouput of fitGeneralizedLotkaVolteraMaximumLikelihood?")
  }

  # Check we can find birth-death parameters
  if ( class(fitted_prey_init) != "numeric" || class(fitted_pred_init) != "numeric" || class(fitted_sigma1) != "numeric" || class(fitted_sigma2) != "numeric") {
    stop("Error in visualizeGeneralizedLVModelFit, cannot find error terms or initial population size parameters in fitted.model, is this argument the ouput of fitGeneralizedLotkaVolteraMaximumLikelihood?")
  }
  
  # Check that the fitted model object contains the right number of parameters
  assertAllFunctionsPassAllTests(preyBirthFxn,
                                 preyDeathFxn,
                                 predBirthFxn,
                                 predDeathFxn,
                                 c(length(fitted_prey_birth),length(fitted_prey_death),length(fitted_pred_birth),length(fitted_pred_death)))
  
  # Simulate trajectories from model
  pred_sim <- matrix(NA,nrow=n_time_points,ncol=n.sim)
  prey_sim <- matrix(NA,nrow=n_time_points,ncol=n.sim)
  pb <- txtProgressBar(min = 0, max = n.sim, style = 3) # report progress
  for (i in 1:n.sim) {
    # Simulate true trajectory
    sim_traj <- sampleGeneralizedLKTrajectory(prey.initial.size=fitted_prey_init,
                                              pred.initial.size=fitted_pred_init,
                                              preyBirthFxn=preyBirthFxn,
                                              preyDeathFxn=preyDeathFxn,
                                              predBirthFxn=predBirthFxn,
                                              predDeathFxn=predDeathFxn,
                                              prey.birth.parameters=fitted_prey_birth,
                                              prey.death.parameters=fitted_prey_death,
                                              pred.birth.parameters=fitted_pred_birth,
                                              pred.death.parameters=fitted_pred_death,
                                              sampled.times=data$time,
                                              n.steps=n_steps)
    # add measurement error
    prey_sim[,i] <- rlnorm(n_time_points,log(sim_traj$prey_trajectory),fitted_sigma1)
    pred_sim[,i] <- rlnorm(n_time_points,log(sim_traj$pred_trajectory),fitted_sigma2)
    setTxtProgressBar(pb, i) # report progress
  }
  
  ## Plot model
  
  # If sigmas and population sizes are too large, user may want to limit y-axis
  # Here we allow specification of ylim as proportion of maximum observed y
  if ( hasArg(relative.y.limit) ) {
    ymax_prey <- list(...)$relative.y.limit*max(data$prey)
    ymax_pred <- list(...)$relative.y.limit*max(data$pred)
  } else {
    ymax_prey <- max(prey_sim)
    ymax_pred <- max(pred_sim)
  }
  
  
  par(mfrow=c(2,1),mai=c(0.8,0.8,0.2,0.2),xpd=TRUE)
  matplot(x=data$time,prey_sim,col=prey_color,type="l",lty=1,xlab="time",ylab="prey pop size",ylim=c(0,ymax_prey))
  lines(data$time,data$prey,lwd=3)
  
  matplot(x=data$time,pred_sim,col=pred_color,type="l",lty=1,xlab="time",ylab="predator pop size",ylim=c(0,ymax_pred))
  lines(data$time,data$pred,lwd=3)
  par(mfrow=c(1,1))
  
}

# This function is called by fitGeneralizedLotkaVolteraMaximumLikelihood and checks the input birth and death rate functions for validity
# Having 80 lines of error checking is cumbersome in the middle of an already long function
assertAllFunctionsPassAllTests <- function(preyBirthFxn,
                                           preyDeathFxn,
                                           predBirthFxn,
                                           predDeathFxn,
                                           n.parameters.in.functions) {  
  # recover()
  # Make sure provided functions are functions
  if ( !class(preyBirthFxn) == "function" ) {
    stop("Argument preyBirthFxn should be the name of a function (and not in quotes).")
  }
  if ( !class(preyDeathFxn) == "function" ) {
    stop("Argument preyDeathFxn should be the name of a function (and not in quotes).")
  }
  if ( !class(predBirthFxn) == "function" ) {
    stop("Argument predBirthFxn should be the name of a function (and not in quotes).")
  }
  if ( !class(predDeathFxn) == "function" ) {
    stop("Argument predDeathFxn should be the name of a function (and not in quotes).")
  }
  
  # Make sure provided functions do not need more than n.parameters.in.functions[i] parameters
  if ( is.na(preyBirthFxn(2,2,rep(2,n.parameters.in.functions[1]))) ) {
    stop("Error with preyBirthFxn. Is n.parameters.in.functions[1] the number of input parameters to this function?")
  }
  if ( is.na(preyDeathFxn(2,2,rep(2,n.parameters.in.functions[2]))) ) {
    stop("Error with preyDeathFxn. Is n.parameters.in.functions[2] the number of input parameters to this function?")
  }
  if ( is.na(predBirthFxn(2,2,rep(2,n.parameters.in.functions[3]))) ) {
    stop("Error with predBirthFxn. Is n.parameters.in.functions[3] the number of input parameters to this function?")
  }
  if ( is.na(predDeathFxn(2,2,rep(2,n.parameters.in.functions[4]))) ) {
    stop("Error with predDeathFxn. Is n.parameters.in.functions[4] the number of input parameters to this function?")
  }
  
  # Try to make sure provided functions are function do not need fewer than n.parameters.in.functions[i] parameters
  # If we have listed the right number of parameters, and we supply one fewer, the result should be an NA
  if ( !is.na(preyBirthFxn(2,2,rep(2,n.parameters.in.functions[1]-1))) ) {
    stop("Error with preyBirthFxn. Is n.parameters.in.functions[1] the number of input parameters to this function?")
  }
  if ( !is.na(preyDeathFxn(2,2,rep(2,n.parameters.in.functions[2]-1))) ) {
    stop("Error with preyDeathFxn. Is n.parameters.in.functions[2] the number of input parameters to this function?")
  }
  if ( !is.na(predBirthFxn(2,2,rep(2,n.parameters.in.functions[3]-1))) ) {
    stop("Error with predBirthFxn. Is n.parameters.in.functions[3] the number of input parameters to this function?")
  }
  if ( !is.na(predDeathFxn(2,2,rep(2,n.parameters.in.functions[4]-1))) ) {
  stop("Error with predDeathFxn. Is n.parameters.in.functions[4] the number of input parameters to this function?")
  }

  # Make sure provided functions return a single numeric
  if ( !is.numeric(preyBirthFxn(2,2,rep(2,n.parameters.in.functions[1]))) || length(preyBirthFxn(2,2,rep(2,n.parameters.in.functions[1]))) != 1 ) {
    stop("Error with preyBirthFxn. Function should return a single numeric value.")
  }
  if ( !is.numeric(preyDeathFxn(2,2,rep(2,n.parameters.in.functions[2])))  || length(preyDeathFxn(2,2,rep(2,n.parameters.in.functions[2]))) != 1 ) {
    stop("Error with preyDeathFxn. Function should return a single numeric value.")
  }
  if ( !is.numeric(predBirthFxn(2,2,rep(2,n.parameters.in.functions[3])))  || length(predBirthFxn(2,2,rep(2,n.parameters.in.functions[3]))) != 1 ) {
    stop("Error with predBirthFxn. Function should return a single numeric value.")
  }
  if ( !is.numeric(predDeathFxn(2,2,rep(2,n.parameters.in.functions[4])))  || length(predDeathFxn(2,2,rep(2,n.parameters.in.functions[4]))) != 1 ) {
    stop("Error with predDeathFxn. Function should return a single numeric value.")
  }
  
  # Make sure functions take arguments with the correct names (this also catches errors that slip through the previous checks)
  if ( class(try(preyBirthFxn(pred.size=1,prey.size=1,prey.birth.parameters=rep(2,n.parameters.in.functions[1])),silent=TRUE)) == "try-error" ) {
    stop("Argument preyBirthFxn. should be a function that takes 2 + n.parameters.in.functions[1] arguments.The arguments must be named pred.size, prey.size, and prey.birth.parameters.")
  }
  if ( class(try(preyDeathFxn(pred.size=1,prey.size=1,prey.death.parameters=rep(2,n.parameters.in.functions[2])),silent=TRUE)) == "try-error" ) {
    stop("Argument preyDeathFxn should be a function that takes 2 + n.parameters.in.functions[2] arguments.The arguments must be named pred.size, prey.size, and prey.death.parameters.")
  }
  if ( class(try(predBirthFxn(pred.size=1,prey.size=1,pred.birth.parameters=rep(2,n.parameters.in.functions[3])),silent=TRUE)) == "try-error" ) {
    stop("Argument predBirthFxn should be a function that takes 2 + n.parameters.in.functions[3] arguments.The arguments must be named pred.size, prey.size, and pred.birth.parameters.")
  }
  if ( class(try(predDeathFxn(pred.size=1,prey.size=1,pred.death.parameters=rep(2,n.parameters.in.functions[4])),silent=TRUE)) == "try-error" ) {
    stop("Argument predDeathFxn should be a function that takes 2 + n.parameters.in.functions[4] arguments.The arguments must be named pred.size, prey.size, and pred.death.parameters.")
  }
  
}
