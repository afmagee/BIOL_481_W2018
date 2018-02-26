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
fit.lv <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=simulated.data,
                                                      preyBirthFxn=preyBirthLV,
                                                      preyDeathFxn=preyDeathLV,
                                                      predBirthFxn=predBirthLV,
                                                      predDeathFxn=predDeathLV,
                                                      n.parameters.in.functions=n.bd.params.lk,
                                                      n.initialization.attempts=1000)

# What are our results for this model?
fit.lv

## Fit a model like Lotka-Voltera, but where the prey birth rate declines as the population gets larger

# This means we have the following numbers of parameters in the birth and death rate functions
n.bd.params.dd <- c(
  2, # Our model now has a rate for small populations and a carrying capacity, 2 parameters
  1, # Lotka-Voltera has a single parameter to describe the prey death rate
  1, # Lotka-Voltera has a single parameter to describe the predator birth rate
  1  # Lotka-Voltera has a single parameter to describe the predator death rate
)

# Now we fit it
fit.dd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=simulated.data,
                                                      preyBirthFxn=preyBirthDensityDependent,
                                                      preyDeathFxn=preyDeathLV,
                                                      predBirthFxn=predBirthLV,
                                                      predDeathFxn=predDeathLV,
                                                      n.parameters.in.functions=n.bd.params.dd,
                                                      n.initialization.attempts=1000)

# What are our results for this model?
fit.dd

## Now we have two competing explanations for our population dynamics. We should compare them
compareModels(list(fit.lv,fit.dd))


## Did either of them do a decent job of recovering the real dynamics?

# Let's look at the first model
visualizeGeneralizedLVModelFit(data=simulated.data,
                               preyBirthFxn=preyBirthLV,
                               preyDeathFxn=preyDeathLV,
                               predBirthFxn=predBirthLV,
                               predDeathFxn=predDeathLV,
                               fitted.model=fit.lv)

# Let's look at the second model
visualizeGeneralizedLVModelFit(data=simulated.data,
                               preyBirthFxn=preyBirthDensityDependent,
                               preyDeathFxn=preyDeathLV,
                               predBirthFxn=predBirthLV,
                               predDeathFxn=predDeathLV,
                               fitted.model=fit.dd)

# Before you hack this script to run your analyses, run this next line. 
# It will clear everything we just did so you don't accidentally mix up simulated and real data
rm(list=ls())