# This function simulates from the Poisson model
# You do have to give it:
#           1) lambda: the "rate" parameter of the Poisson, not interpretable here as a mutation rate, but as the product of the mutation rate and the number of plated bacteria
#           2) real.data: your data
#           3) num.draws: number of simulated experiments (each of equal size to the real one)
# The summary statistics are calculated for each simulated experiment
#           1) mean: the mean number of surviving bacteria
#           2) variance: the variance deviation of the number of surviving bacteria 
#           3) relative.variance: the variance of the number of surviving bacteria divided by the mean (should be 1 for a Poisson model)
#           4) kurtosis: a measure of the tendency of the model to display particularly large values (numbers of surviving colonies)

poissonModel <- function(lambda, real.data, num.draws=10000) {
  
  # recover()
  
  # Draw samples from the specified Poisson distribution
  # We're exploiting the fact that each datapoint is iid (independent and identically distributed)
  # This means we can simulate all the data points, and arbitrarily assign them to the experiment they came from
  poisson_draws <- rpois(num.draws*length(real.data),lambda)
  
  # Assemble these into a number of replicate datasets
  y_rep <- matrix(poisson_draws,nrow=length(real.data),ncol=num.draws,byrow=TRUE)
  
  simulationVisualSummary(y_rep, real.data)

}

# This function does a statistical analysis on your behalf!
# Thus, you do not need to give it a value to control the Poisson distribution, it finds that value.
# Everything else works as in poissonModel()
poissonBayes <- function(real.data, num.draws=10000, summary.statistic=c("mean","standard.deviation","relative.variance","mean.absolute.deviation")[1], ...) {
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
  posterior_gamma_shape <- gamma_prior_shape + sum(real.data)
  posterior_gamma_rate  <- gamma_prior_rate + length(real.data)
  
  # Draw from posterior of Poisson rate parameter
  poisson_rate_draws <- rgamma(num.draws, shape=posterior_gamma_shape, rate=posterior_gamma_rate)
  
  # Simulate datasets according to this posterior ("draw replicate datasets")
  # We're exploiting the fact that each datapoint is iid (independent and identically distributed)
  # This means we can simulate all the data points, and arbitrarily assign them to the experiment they came from
  y_rep <- rpois(num.draws*length(real.data),poisson_rate_draws)
  
  # Fold replicates into a table
  y_rep <- matrix(y_rep,nrow=length(real.data),ncol=num.draws,byrow=TRUE)
  
  simulationVisualSummary(y_rep, real.data)
}

# simulated.data is the simulated data, as a matrix (if you call this on any data generated here, it will be appropriately formatted)
# real.data is the real data as a vector
simulationVisualSummary <- function(simulated.data,real.data) {

  # recover()
  
  # Compute all summary statistics
  sim_means <- colMeans(simulated.data)
  observed_mean <- mean(real.data)
  
  sim_vars <- apply(simulated.data,2,var)
  observed_var <- var(real.data)
  
  sim_rel_vars <- apply(simulated.data,2,var) / colMeans(simulated.data)
  observed_rel_var <- var(real.data) / mean(real.data)
  
  # Skewness is not as informative here as I imagined
  # sim_nonparametric_skews <- colMeans(simulated.data) - apply(simulated.data,2,median)
  # observed_nonparametric_skew <- mean(simulated.data) - median(simulated.data)
  
  sim_kurtosis <- apply(simulated.data,2,kurtosis)
  observed_kurtosis <- kurtosis(real.data)
  
  # Calculate the posterior-predictive p-values
  p_mean <- sum(sim_means > observed_mean)/length(sim_means)
  p_var <- sum(sim_vars > observed_var)/length(sim_vars)
  p_rel_var <- sum(sim_rel_vars > observed_rel_var)/length(sim_rel_vars)
  # p_skew <- sum(sim_nonparametric_skews > observed_nonparametric_skew)/length(sim_nonparametric_skews)
  p_kurt <- sum(sim_kurtosis > observed_kurtosis)/length(sim_kurtosis)
  cat("Probability of seeing value of mean              greater than observed: ",p_mean,"\n",sep="")
  cat("Probability of seeing value of variance          greater than observed: ",p_var,"\n",sep="")
  cat("Probability of seeing value of relative variance greater than observed: ",p_rel_var,"\n",sep="")
  # cat("Probability of seeing value of skewness          greater than observed: ",p_skew,"\n",sep="")
  cat("Probability of seeing value of kurtosis          greater than observed: ",p_kurt,"\n",sep="")
  
  # Prepare for plotting, is our observed value even in the range of simulated values?
  obs_mean_in_range <- min(sim_means) < observed_mean && max(sim_means) > observed_mean
  obs_var_in_range <- min(sim_vars) < observed_var && max(sim_vars) > observed_var
  obs_rel_var_in_range <- min(sim_rel_vars) < observed_rel_var && max(sim_rel_vars) > observed_rel_var
  # obs_nonparametric_skew_in_range <- min(sim_nonparametric_skews) < observed_nonparametric_skew && max(sim_nonparametric_skews) > observed_nonparametric_skew
  obs_kurtosis_in_range <- min(sim_kurtosis) < observed_kurtosis && max(sim_kurtosis) > observed_kurtosis
  
  # Start plot, save as an object (weird, right?)
  plot_mean <- ggplot() + 
    # Tell it about where the data is
    aes(sim_means) + 
    # Make sure we can see the real value, even if it's far away
    expand_limits(x=observed_mean) + 
    # Titles and axis labels
    ggtitle("Distribution of test statistic") + xlab("mean") + 
    # The actual histogram (with on-the-fly number of bins)
    geom_histogram(breaks=pretty(sim_means,n=ifelse(obs_mean_in_range,50,500),min.n=ifelse(obs_mean_in_range,40,400)),fill="deepskyblue2") +
    # The real value
    geom_vline(xintercept=observed_mean,color="darkorange2",lwd=1.2)
  
  plot_var <- ggplot() + 
    aes(sim_vars) + 
    expand_limits(x=observed_var) + 
    ggtitle("Distribution of test statistic") + 
    xlab("variance") + 
    geom_histogram(bins=ifelse(obs_var_in_range,50,500),fill="deepskyblue2") +
    geom_vline(xintercept=observed_var,color="darkorange2",lwd=1.2)
  
  plot_rel_var <- ggplot() + 
    aes(sim_rel_vars) + 
    expand_limits(x=observed_rel_var) + 
    ggtitle("Distribution of test statistic") + 
    xlab("variance/mean") + 
    geom_histogram(bins=ifelse(obs_rel_var_in_range,50,500),fill="deepskyblue2") +
    geom_vline(xintercept=observed_rel_var,color="darkorange2",lwd=1.2)
  
  # plot_nonparametric_skew <- ggplot() + 
  #   aes(sim_nonparametric_skews) + 
  #   expand_limits(x=observed_nonparametric_skew) + 
  #   ggtitle("Distribution of test statistic") + 
  #   xlab("skewness") + 
  #   geom_histogram(bins=ifelse(obs_nonparametric_skew_in_range,50,500),fill="deepskyblue2") +
  #   geom_vline(xintercept=observed_nonparametric_skew,color="darkorange2",lwd=1.2)

  plot_kurtosis <- ggplot() + 
    aes(sim_kurtosis) + 
    expand_limits(x=observed_kurtosis) + 
    ggtitle("Distribution of test statistic") + 
    xlab("kurtosis") + 
    geom_histogram(bins=ifelse(obs_kurtosis_in_range,25,250),fill="deepskyblue2") +
    geom_vline(xintercept=observed_kurtosis,color="darkorange2",lwd=1.2)
  
  grid.arrange(plot_mean, plot_var, plot_rel_var, plot_kurtosis, ncol=2)
}

# # Test on some data from the original Luria-Delbruck experiment!
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"mean")
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"standard.deviation")
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"relative.variance")
# poissonBayes(c(1,1,3,5,5,6,35,64,107,rep(0,11)),1000,"mean.absolute.deviation")

kurtosis <- function(x) {
  mean((x - mean(x))^4)/(mean((x - mean(x))^2)^2)
}
