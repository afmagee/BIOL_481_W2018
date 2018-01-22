# This function simulates from the Poisson model
# You do have to give it:
#           1) lambda: the "rate" parameter of the Poisson, not interpretable here as a mutation rate, but as the product of the mutation rate and the number of plated bacteria
#           2) data.points: your data
#           3) num.draws: number of simulated experiments (each of equal size to the real one)
#           4) summary.statistic: how should the data (for each simulated experiment) be summarized
# The summary statistics are calculated for each simulated experiment
#           1) mean: the mean number of surviving bacteria
#           2) standard.deviation: the standard deviation of the number of surviving bacteria 
#           3) relative.variance: the variance of the number of surviving bacteria divided by the mean (should be 1 for a Poisson model)
#           4) mean.absolute.deviation: like the standard deviation or variance, but measuring absolute differences from the mean instead of squared distances

poissonModel <- function(lambda, data.points, num.draws=10000, summary.statistic=c("mean","standard.deviation","relative.variance","mean.absolute.deviation")[1]) {

  # Draw samples from the specified Poisson distribution
  # We're exploiting the fact that each datapoint is iid (independent and identically distributed)
  # This means we can simulate all the data points, and arbitrarily assign them to the experiment they came from
  poisson_draws <- rpois(num.draws*length(data.points),lambda)
  
  # Assemble these into a number of replicate datasets
  y_rep <- matrix(poisson_draws,nrow=length(data.points),ncol=num.draws,byrow=TRUE)
  
  # Compute desired summary statistic
  if ( summary.statistic == "mean" ) {
    summary_stats <- colMeans(y_rep)
    observed_stat <- mean(data.points)
  } else if ( summary.statistic == "standard.deviation" ) {
    summary_stats <- apply(y_rep,2,sd)
    observed_stat <- sd(data.points)
  } else if ( summary.statistic == "relative.variance" ) {
    summary_stats <- apply(y_rep,2,var) / colMeans(y_rep)
    observed_stat <- var(data.points) / mean(data.points)
  } else if ( summary.statistic == "mean.absolute.deviation" ) {
    summary_stats <- colMeans(abs(y_rep-colMeans(y_rep)))
    observed_stat <- mean(abs(data.points - mean(data.points)))
  } else {
    stop("Option for argumant \"summary.statistic\" not recognized")
  }
  
  # Calculate the posterior-predictive p-value
  p <- sum(summary_stats > observed_stat)/num.draws
  cat("Probability of seeing value of ",summary.statistic," greater than observed: ",p,"\n",sep="")

  # Pre-plotting computation, is our observed value even in the range of simulated values?
  obs_in_range <- min(summary_stats) < observed_stat && observed_stat < max(summary_stats)
  
  # Start plot
  ggplot() + # We can break a line into multiple lines at math operators
    
    # Tell it about where the data is
    aes(summary_stats) + 
    
    # Avoid 
    expand_limits(x=observed_stat) + 
    # 
    
    # Titles and axis labels
    ggtitle("Distribution of test statistic") + xlab(summary.statistic) + 
    
    # Make a decision on the number of bins
    # If our observed value is within the simulated range, we'll use 50 bins
    # If not, we need way more bins to get any resolution on the distribution
    geom_histogram(bins=ifelse(obs_in_range,50,500),fill="deepskyblue2") +
    
    # Add a line
    geom_vline(xintercept=observed_stat,color="darkorange2",lwd=1.2)
  
}

# This function does a statistical analysis on your behalf!
# Thus, you do not need to give it a value to control the Poisson distribution, it finds that value.
# Everything else works as in poissonModel()
poissonBayes <- function(data.points, num.draws=10000, summary.statistic=c("mean","standard.deviation","relative.variance","mean.absolute.deviation")[1], ...) {
  # Check for user input of prior (gamma) distribution
  # The prior can be thought of as follows: I have observed an additional gamma_prior_rate data points with a sum of gamma_prior_shape
  # This can throw estimation off if gamma_prior_rate is large compared to the real number of observations and the real mean is large
  if (hasArg(gamma.shape)) {
    gamma_prior_shape <- list(...)$gamma.shape
  } else {
    gamma_prior_shape <- 1e-10
  }
  if (hasArg(gamma.rate)) {
    gamma_prior_rate <- list(...)$gamma.rate
  } else {
    gamma_prior_rate <- 1e-10
  }
  
  # Compute the posterior distribution
  posterior_gamma_shape <- gamma_prior_shape + sum(data.points)
  posterior_gamma_rate  <- gamma_prior_rate + length(data.points)
  
  # Draw from posterior of Poisson rate parameter
  poisson_rate_draws <- rgamma(num.draws, shape=posterior_gamma_shape, rate=posterior_gamma_rate)
  
  # Simulate datasets according to this posterior ("draw replicate datasets")
  # We're exploiting the fact that each datapoint is iid (independent and identically distributed)
  # This means we can simulate all the data points, and arbitrarily assign them to the experiment they came from
  y_rep <- rpois(num.draws*length(data.points),poisson_rate_draws)
  
  # Fold replicates into a table
  y_rep <- matrix(y_rep,nrow=length(data.points),ncol=num.draws,byrow=TRUE)
  
  # recover()
  
  # Compute desired summary statistic
  if ( summary.statistic == "mean" ) {
    summary_stats <- colMeans(y_rep)
    observed_stat <- mean(data.points)
  } else if ( summary.statistic == "standard.deviation" ) {
    summary_stats <- apply(y_rep,2,sd)
    observed_stat <- sd(data.points)
  } else if ( summary.statistic == "relative.variance" ) {
    summary_stats <- apply(y_rep,2,var) / colMeans(y_rep)
    observed_stat <- var(data.points) / mean(data.points)
  } else if ( summary.statistic == "mean.absolute.deviation" ) {
    summary_stats <- colMeans(abs(y_rep-colMeans(y_rep)))
    observed_stat <- mean(abs(data.points - mean(data.points)))
  } else {
    stop("Option for argumant \"summary.statistic\" not recognized")
  } 
  # recover()
  
  # Calculate the posterior-predictive p-value
  pppval <- sum(summary_stats > observed_stat)/num.draws
  cat("Posterior-predictive p-value:",pppval,"\n")
  
  # Pre-plotting computation, is our observed value even in the range of simulated values?
  obs_in_range <- min(summary_stats) < observed_stat && observed_stat < max(summary_stats)
  
  # Start plot
  ggplot() + # We can break a line into multiple lines at math operators
    
    # Tell it about where the data is
    aes(summary_stats) + 
    
    # Avoid 
    expand_limits(x=observed_stat) + 
    # 
    
    # Titles and axis labels
    ggtitle("Posterior-predictive distribution") + xlab(summary.statistic) + 
    
    # Make a decision on the number of bins
    # If our observed value is within the simulated range, we'll use 50 bins
    # If not, we need way more bins to get any resolution on the distribution
    geom_histogram(bins=ifelse(obs_in_range,50,500),fill="deepskyblue2") +
    
    # Add a line
    geom_vline(xintercept=observed_stat,color="darkorange2",lwd=1.2)
  
}

# # Test on some data from the original Luria-Delbruck experiment!
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"mean")
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"standard.deviation")
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"relative.variance")
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"mean.absolute.deviation")
