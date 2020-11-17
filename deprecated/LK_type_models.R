###################################################
# Models where predators and prey are independent #
###################################################

independentExponentialGrowth <- function(prey.size,pred.size,theta) {
  # This is standard exponential growth of both populations individually
  # We label the parameters as follows:
  #      theta[1] = (birth - death) rate of the prey (not inough information in this world to estimate birth and death rates separately here)
  #      theta[2] = (birth - death) rate of the predator (not inough information in this world to estimate birth and death rates separately here)
  prey_slope <- theta[1]*prey.size
  pred_slope <- theta[2]*pred.size
  slope_vector <- c(prey_slope,pred_slope)
  return(slope_vector)
}

independentLogisticGrowth <- function(prey.size,pred.size,theta) {
  # This is standard exponential growth of both populations individually
  # We label the parameters as follows:
  #      theta[1] = (birth - death) rate of the prey (not inough information in this world to estimate birth and death rates separately here)
  #      theta[2] = carrying capacity of prey
  #      theta[3] = (birth - death) rate of the predator (not inough information in this world to estimate birth and death rates separately here)
  #      theta[4] = carrying capacity of predator
  prey_slope <- theta[1]*(1 - prey.size/theta[2])
  pred_slope <- theta[3]*(1 - pred.size/theta[4])
  slope_vector <- c(prey_slope,pred_slope)
  return(slope_vector)
}

###############################################################################
# Lotka-Voltera and similar models without density-dependent population sizes #
###############################################################################

standardLK <- function(prey.size,pred.size,theta) {
  # This is just the standard Lotka-Voltera machinary
  # We label the parameters as follows:
  #      theta[1] = birth rate of the prey 
  #      theta[2] = death rate of the prey (depends on predator population size)
  #      theta[3] = death rate of the predator 
  #      theta[4] = birth rate of the predator (depends on prey population size)
  prey_slope <- (theta[1] - theta[2]*pred.size)*prey.size
  pred_slope <- (-theta[3] + theta[4]*prey.size)*pred.size
  slope_vector <- c(prey_slope,pred_slope)
  return(slope_vector)
}

constantAttackLK <- function(prey.size,pred.size,theta) {
  # This is Lotka-Voltera-like machinary with what May (1976) calls a constant attack rate formulation
  # We label the parameters as follows:
  #      theta[1] = birth rate of the prey 
  #      theta[2] = death rate of the prey (here does not depend on prey population size)
  #      theta[3] = death rate of the predator 
  #      theta[4] = birth rate of the predator (depends on prey population size)
  prey_slope <- theta[1]*prey.size - theta[2]*pred.size
  pred_slope <- (-theta[3] + theta[4]*prey.size)*pred.size
  slope_vector <- c(prey_slope,pred_slope)
  return(slope_vector)
}

ivlevLK <- function(prey.size,pred.size,theta) {
  # This is Lotka-Voltera-like machinary with what May (1976) calls a Holling Type II (invertebrate) response due to Ivlev
  # We label the parameters as follows:
  #      theta[1] = birth rate of the prey 
  #      theta[2] = death rate of the prey when predator population size is very large
  #      theta[3] = describes how the death rate of the prey depends on the population size (a density-dependent effect)
  #      theta[4] = death rate of the predator 
  #      theta[5] = birth rate of the predator (depends on prey population size)
  prey_slope <- (theta[1] - theta[2]*(1 - exp(-theta[3]*pred.size)))*prey.size
  pred_slope <- (-theta[4] + theta[5]*prey.size)*pred.size
  slope_vector <- c(prey_slope,pred_slope)
  return(slope_vector)
}

#####################################################################
# Lotka-Voltera-like models with density-dependent population sizes #
#####################################################################

predDensityRegulatedLK <- function(prey.size,pred.size,theta) {
  # This is Lotka-Voltera-like machinary with logistic growth of the predator population size from May (1976)
  # We label the parameters as follows:
  #      theta[1] = birth rate of the prey 
  #      theta[2] = death rate of the prey (depends on predator population size)
  #      theta[3] = death rate of the predator 
  #      theta[4] = birth rate of the predator when predator population is small
  #      theta[5] = describes how the predator carrying capacity relates to the prey population size
  prey_slope <- (theta[1] - theta[2]*pred.size)*prey.size
  pred_slope <- (-theta[3] + theta[4]*(1 - theta[5]*pred.size/prey.size))*pred.size
  slope_vector <- c(prey_slope,pred_slope)
}

preyDensityRegulatedLK <- function(prey.size,pred.size,theta) {
  # This is Lotka-Voltera-like machinary with logistic growth of the predator population size from May (1976)
  # We label the parameters as follows:
  #      theta[1] = birth rate of the prey when prey population size is small
  #      theta[2] = prey population carrying capacity
  #      theta[3] = death rate of the prey (depends on predator population size)
  #      theta[4] = death rate of the predator 
  #      theta[5] = birth rate of the predator (depends on prey population size)
  prey_slope <- (theta[1]*(1 - prey.size/theta[2]) - theta[3]*pred.size)*prey.size
  pred_slope <- (-theta[4] + theta[5]*prey.size)*pred.size
  slope_vector <- c(prey_slope,pred_slope)
}

predAndPreyDensityRegulatedLK <- function(prey.size,pred.size,theta) {
  # This is Lotka-Voltera-like machinary with logistic growth of the predator population size from May (1976)
  # We label the parameters as follows:
  #      theta[1] = birth rate of the prey when prey population size is small
  #      theta[2] = prey population carrying capacity
  #      theta[3] = death rate of the prey (depends on predator population size)
  #      theta[4] = death rate of the predator 
  #      theta[5] = birth rate of the predator when predator population is small
  #      theta[6] = describes how the predator carrying capacity relates to the prey population size
  prey_slope <- (theta[1]*(1 - prey.size/theta[2]) - theta[3]*pred.size)*prey.size
  pred_slope <- (-theta[4] + theta[5]*(1 - theta[6]*pred.size/prey.size))*pred.size
  slope_vector <- c(prey_slope,pred_slope)
}



# Figure out how to model these things with birth and death rates
# Put all the pop sizes on the same scale
# In intro/worksheet, explain the coefficients