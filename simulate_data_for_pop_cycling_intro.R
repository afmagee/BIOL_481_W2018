set.seed(8472)

simulated.data <- sampleGeneralizedLKTrajectory(
  40,
  2,
  preyBirthLK,
  preyDeathLK,
  predBirthLK,
  predDeathLK,
  c(0.9),
  c(0.01),
  c(0.01),
  c(0.9),
  0:20,
  10000
)

simulated.data$times <- c(1:21)

plotPredatorPrey(simulated.data)

sigma.prey <- 0.1
sigma.pred <- 0.1

simulated.data$prey_trajectory <- rlnorm(21,log(simulated.data$prey_trajectory),sigma.prey)
simulated.data$pred_trajectory <- rlnorm(21,log(simulated.data$pred_trajectory),sigma.pred)

plotPredatorPrey(simulated.data)

write.csv(simulated.data,"simulated_data.csv",quote=FALSE,row.names=FALSE)
