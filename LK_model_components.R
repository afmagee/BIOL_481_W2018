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
  # This has Lotka-Voltera-like dynamics if the prey size is very large
  return(prey.death.parameters[1]*(1 - exp(-prey.death.parameters[2]*prey.size)))
}

preyDeathWatt <- function(pred.size,prey.size,prey.death.parameters) {
  # From Chapter 5, May (1976), also called a Holling Type II (Invertebrate) response
  # This is pretty different from Lotka-Voltera, unless the prey size is large and the predator size is small
  return(prey.death.parameters[1]*(1 - exp(-prey.death.parameters[2]*prey.size*(pred.size^(1-prey.death.parameters[3])))))
}

## Density-dependent population birth rates
predBirthDensityDependent1 <- function(pred.size,prey.size,pred.birth.parameters) {
  # From Chapter 5, May (1976)
  # The predators have a growth rate all their own, but the carrying capacity depends on the prey population size
  return(pred.birth.parameters[1]*(1 - pred.birth.parameters[2]*pred.size/prey.size))
}

predBirthDensityDependent2 <- function(pred.size,prey.size,pred.birth.parameters) {
  # This has Lotka-Voltera dynamics when predator population size is small
  # Growth rate slows as population increases
  return(prey.size*pred.birth.parameters[1]*(1 - pred.size/pred.birth.parameters[2]))
}

preyBirthDensityDependent <- function(pred.size,prey.size,prey.birth.parameters) {
  # Growth rate slows as population increases
  return(prey.birth.parameters[1]*(1 - prey.size/prey.birth.parameters[2]))
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


system.time({lh.dddd <- 
  fitGeneralizedLotkaVolteraMaximumLikelihood(lh.data,
                                            preyBirthDensityDependent,
                                            preyDeathLK,
                                            predBirthDensityDependent2,
                                            predDeathLK,
                                            c(2,1,2,1))
})

mm.data <- read.csv("population_cycling/mink_muskrat.csv")

system.time({mm.wapreddd <- 
  fitGeneralizedLotkaVolteraMaximumLikelihood(mm.data,
                                              preyBirthLK,
                                              preyDeathWatt,
                                              predBirthDensityDependent1,
                                              predDeathLK,
                                              c(1,3,2,1))
})

assessGeneralizedLKModelFit(mm.data,
                            preyBirthLK,
                            preyDeathIvlev,
                            predBirthHolling,
                            predDeathLK,
                            mm.ivho)


fitGeneralizedLotkaVolteraMaximumLikelihood(lh.data,
                                            preyBirthDensityDependent,
                                            preyDeathLK,
                                            predBirthLK,
                                            predDeathLK,
                                            c(1,1,1,1))
