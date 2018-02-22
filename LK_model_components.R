# ALL OF THESE FUNCTIONS ARE PER-CAPITA (because that is the way they will be interpreted in inference)

## The components of the basic Lotka-Voltera model
predBirthLK <- function(pred.size,prey.size,pred.birth.parameters) {
  return(prey.size*pred.birth.parameters[1])
}

predDeathLK <- function(pred.size,prey.size,pred.death.parameters) {
  return(pred.death.parameters[1])
}

preyBirthLK <- function(pred.size,prey.size,prey.birth.parameters) {
  return(prey.birth.parameters[1])
}

preyDeathLK <- function(pred.size,prey.size,prey.death.parameters) {
  return(pred.size*prey.death.parameters[1])
}

## Different takes on prey death
preyDeathIvlev <- function(pred.size,prey.size,prey.death.parameters) {
  # From Chapter 5, May (1976), also called a Holling Type II (Invertebrate) response
  # 
  return(prey.death.parameters[1]*(1 - exp(-prey.death.parameters[2]*pred.size)))
}


## Density-dependent population sizes
predBirthDensityDependent <- function(pred.size,prey.size,pred.birth.parameters) {
  # From Chapter 5, May (1976), also called a Holling Type II (Invertebrate) response
  return(pred.birth.parameters[1]*(1 - pred.birth.parameters[2]*pred.size/prey.size))
}

preyBirthDensityDependent <- function(pred.size,prey.size,prey.birth.parameters) {
  return(prey.birth.parameters[1]*(1 - prey.size/prey.birth.parameters[2]))
}



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


simdata <- 
sampleGeneralizedLKTrajectory(34,6,
                              preyBirthLK,
                              preyDeathLK,
                              predBirthLK,
                              predDeathLK,
                              prey.birth.parameters=c(0.5),
                              prey.death.parameters=c(0.3),
                              pred.birth.parameters=c(0.2),
                              pred.death.parameters=c(0.8),
                              0:10,
                              5000
                              )
simdata$time <- 0:10

preyBirthLK(34,6,0.5) - preyDeathLK(34,6,0.3)

logLikelihoodGeneralizedLK(data=lh.data,
                           42,
                           5,
                           preyBirthLK,
                           preyDeathLK,
                           predBirthLK,
                           predDeathLK,
                           c(1),
                           c(1),
                           c(1),
                           c(1),
                           0.25,
                           0.25,
                           1000,
                           FALSE)

fitGeneralizedLotkaVolteraMaximumLikelihood(lh.data,
                                            preyBirthLK,
                                            preyDeathLK,
                                            predBirthLK,
                                            predDeathLK,
                                            c(1,1,1,1))


system.time({lh.iv <- 
  fitGeneralizedLotkaVolteraMaximumLikelihood(lh.data,
                                            preyBirthLK,
                                            preyDeathIvlev,
                                            predBirthLK,
                                            predDeathLK,
                                            c(1,2,1,1))
})



fitGeneralizedLotkaVolteraMaximumLikelihood(lh.data,
                                            preyBirthDensityDependent,
                                            preyDeathLK,
                                            predBirthLK,
                                            predDeathLK,
                                            c(1,1,1,1))
