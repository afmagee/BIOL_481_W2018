## Source the functions we need
source("population_inference.R")
source("LK_model_components.R")

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
fit.lk <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=simulated.data,
                                                      preyBirthFxn=preyBirthLK,
                                                      preyDeathFxn=preyDeathLK,
                                                      predBirthFxn=predBirthLK,
                                                      predDeathFxn=predDeathLK,
                                                      n.parameters.in.functions=n.bd.params.lk)
})
