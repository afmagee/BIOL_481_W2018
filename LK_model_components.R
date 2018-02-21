# ALL OF THESE FUNCTIONS ARE PER-CAPITA (because that is the way they will be interpreted in inference)

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


sampleGeneralizedLKTrajectory(40,4,
                              preyBirthLK,
                              preyDeathLK,
                              predBirthLK,
                              predDeathLK,
                              c(1),
                              c(1),
                              c(1),
                              c(1),
                              0:10,
                              1000
                              )