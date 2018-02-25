## Source the functions we need
source("population_inference.R")
source("LV_model_components.R")

## Read in the data
simulated.data <- read.csv("simulated_data.csv")

## Plot the observed trajectories
plotPredatorPrey(simulated.data)

## Fit a standard Lotka-Voltera model

# This means we have the following numbers of parameters in the birth and death rate functions
n.bd.params.lk <- c(
  1, # Lotka-Voltera has a single parameter to describe the prey birth rate
  1, # Lotka-Voltera has a single parameter to describe the prey death rate
  1, # Lotka-Voltera has a single parameter to describe the predator birth rate
  1  # Lotka-Voltera has a single parameter to describe the predator death rate
)

# Now we fit it
system.time({
fit.lv <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=simulated.data,
                                                      preyBirthFxn=preyBirthLV,
                                                      preyDeathFxn=preyDeathLV,
                                                      predBirthFxn=predBirthLV,
                                                      predDeathFxn=predDeathLV,
                                                      n.parameters.in.functions=n.bd.params.lk,
                                                      n.initialization.attempts=1000)
})

## Fit a model like Lotka-Voltera, but where the prey birth rate declines as the population gets larger

# This means we have the following numbers of parameters in the birth and death rate functions
n.bd.params.dd <- c(
  2, # Our model now has a rate for small populations and a carrying capacity, 2 parameters
  1, # Lotka-Voltera has a single parameter to describe the prey death rate
  1, # Lotka-Voltera has a single parameter to describe the predator birth rate
  1  # Lotka-Voltera has a single parameter to describe the predator death rate
)

# Now we fit it
system.time({
  fit.dd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=simulated.data,
                                                        preyBirthFxn=preyBirthDensityDependent,
                                                        preyDeathFxn=preyDeathLV,
                                                        predBirthFxn=predBirthLV,
                                                        predDeathFxn=predDeathLV,
                                                        n.parameters.in.functions=n.bd.params.dd,
                                                        n.initialization.attempts=1000)
})

## Now we have two competing explanations for our population dynamics. We should compare them
compareModels(list(fit.lv,fit.dd))


## Did either of them do a decent job of recovering the real dynamics?

# Let's look at the best model
visualizeGeneralizedLVModelFit(simulated.data,
                               preyBirthFxn=preyBirthLV,
                               preyDeathFxn=preyDeathLV,
                               predBirthFxn=predBirthLV,
                               predDeathFxn=predDeathLV,
                               fit.dd)

# Let's look at the next-best model
visualizeGeneralizedLVModelFit(simulated.data,
                               preyBirthFxn=preyBirthDensityDependent,
                               preyDeathFxn=preyDeathLV,
                               predBirthFxn=predBirthLV,
                               predDeathFxn=predDeathLV,
                               fit.lv)

