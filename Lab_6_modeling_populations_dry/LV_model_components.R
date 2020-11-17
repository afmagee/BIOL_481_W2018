# All of these functions are per-capita rates (because that is the way they will be interpreted in inference)
# All of these functions take in the prey and predator sizes, a parameter (or multiple parameters), and return the birth or death rate for those values
# The parameters are named to match what they are for (e.g. pred.birth.parameters for predator birth functions)

## The components of the basic Lotka-Voltera model
predBirthLV <- function(pred.size,prey.size,pred.birth.parameters) {
  # The birth rate of the predator is a linear function of the number of prey
  return(prey.size*pred.birth.parameters[1])
}

predDeathLV <- function(pred.size,prey.size,pred.death.parameters) {
  # The death rate of the predators is constant
  return(pred.death.parameters[1])
}

preyBirthLV <- function(pred.size,prey.size,prey.birth.parameters) {
  # The birth rate of the prey is constant
  return(prey.birth.parameters[1])
}

preyDeathLV <- function(pred.size,prey.size,prey.death.parameters) {
  # The death rate of the prey is a linear function of the number of predators
  return(pred.size*prey.death.parameters[1])
}

## Different models of prey birth rates (here density-dependent)
preyBirthDensityDependent <- function(pred.size,prey.size,prey.birth.parameters) {
  # Growth rate slows as population increases
  return(prey.birth.parameters[1]*(1 - prey.size/prey.birth.parameters[2]))
}

## Different models of prey death rates
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

## Different models of predator birth rates (here density-dependent)
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

