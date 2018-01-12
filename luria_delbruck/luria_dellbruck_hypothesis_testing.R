poissonModel <- function(data.points, lambda, num.draws, summary.statistic=c("mean","variance","coefficient.of.variation","mean.absolute.deviation")[1]) {
  
  # Draw samples from the specified Poisson distribution
  poisson_draws <- rpois(num.draws*length(data.points),lambda)
  
  # Assemble these into a number of replicate datasets
  y_rep <- matrix(poisson_draws,nrow=length(data.points),ncol=num.draws,byrow=TRUE)
  
  # Compute desired summary statistic
  if ( summary.statistic == "mean" ) {
    summary_stats <- colMeans(y_rep)
    observed_stat <- mean(data.points)
  } else if ( summary.statistic == "variance" ) {
    summary_stats <- apply(y_rep,2,var)
    observed_stat <- var(data.points)
  } else if ( summary.statistic == "coefficient.of.variation" ) {
    summary_stats <- apply(y_rep,2,var) / colMeans(y_rep)
    observed_stat <- var(data.points) / mean(data.points)
  } else if ( summary.statistic == "mean.absolute.deviation" ) {
    summary_stats <- colMeans(y_rep-colMeans(y_rep))
    observed_stat <- mean(abs(data.points - mean(data.points)))
  } else {
    stop("Option for argumant \"summary.statistic\" not recognized")
  }
  
  # Plot
  hist(summary_stats,main=paste0("Distribution of summary statistic from Poisson with lambda=",lambda),xlab=summary.statistic,freq=FALSE,xlim=c(min(min(summary_stats),observed_stat*0.9),max((summary_stats),observed_stat*1.1)))
  abline(v=observed_stat,col="red")
  
  # Calculate the posterior-predictive p-value
  p <- sum(summary_stats > observed_stat)/num.draws
  cat("Probability of observing value of",summary.statistic," greater than observed:",p,"\n",sep="")
}
poissonBayes <- function(data.points, num.draws, summary.statistic=c("mean","variance","coefficient.of.variation","mean.absolute.deviation")[1], ...) {
  # Check for user input of prior (gamma) distribution
  if (hasArg(gamma.shape)) {
    gamma_prior_shape <- list(...)$gamma.shape
  } else {
    gamma_prior_shape <- 0.5
  }
  if (hasArg(gamma.rate)) {
    gamma_prior_rate <- list(...)$gamma.rate
  } else {
    gamma_prior_rate <- 0.5
  }
  
  # Compute the posterior distribution
  posterior_gamma_shape <- gamma_prior_shape + sum(data.points)
  posterior_gamma_rate  <- gamma_prior_rate + length(data.points)
  
  # Draw from posterior of Poisson rate parameter
  poisson_rate_draws <- rgamma(num.draws, shape=posterior_gamma_shape, rate=posterior_gamma_rate)
  
  # Simulate datasets according to this posterior ("draw replicate datasets")
  y_rep <- rpois(num.draws*length(data.points),poisson_rate_draws)
  
  # Fold replicates into a table
  y_rep <- matrix(y_rep,nrow=length(data.points),ncol=num.draws,byrow=TRUE)
  
  # recover()
  
  # Compute desired summary statistic
  if ( summary.statistic == "mean" ) {
    summary_stats <- colMeans(y_rep)
    observed_stat <- mean(data.points)
  } else if ( summary.statistic == "variance" ) {
    summary_stats <- apply(y_rep,2,var)
    observed_stat <- var(data.points)
  } else if ( summary.statistic == "coefficient.of.variation" ) {
    summary_stats <- apply(y_rep,2,var) / colMeans(y_rep)
    observed_stat <- var(data.points) / mean(data.points)
  } else if ( summary.statistic == "mean.absolute.deviation" ) {
    summary_stats <- colMeans(y_rep-colMeans(y_rep))
    observed_stat <- mean(abs(data.points - mean(data.points)))
  } else {
    stop("Option for argumant \"summary.statistic\" not recognized")
  } 
  
  # Plot
  hist(summary_stats,main="Posterior-predictive distribution",xlab=summary.statistic,freq=FALSE,xlim=c(min(min(summary_stats),observed_stat*0.9),max((summary_stats),observed_stat*1.1)))
  abline(v=observed_stat,col="red")
  
  # Calculate the posterior-predictive p-value
  pppval <- sum(summary_stats > observed_stat)/num.draws
  cat("Posterior-predictive p-value:",pppval,"\n")
}


poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"mean")
poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"variance")
poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"coefficient.of.variation")
poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"mean.absolute.deviation")

