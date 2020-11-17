simulationVisualSummary <- function(simulated.data,real.data) {
  
  # recover()
  
  # Compute all summary statistics
  sim_means <- colMeans(simulated.data)
  observed_mean <- mean(real.data)
  
  sim_vars <- apply(simulated.data,2,var)
  observed_var <- var(real.data)
  
  sim_rel_vars <- apply(simulated.data,2,var) / colMeans(simulated.data)
  observed_rel_var <- var(real.data) / mean(real.data)
  
  sim_kurtosis <- apply(simulated.data,2,kurtosis)
  observed_kurtosis <- kurtosis(real.data)
  
  # Calculate the posterior-predictive p-values
  p_mean <- sum(sim_means > observed_mean)/length(sim_means)
  p_var <- sum(sim_vars > observed_var)/length(sim_vars)
  p_rel_var <- sum(sim_rel_vars > observed_rel_var)/length(sim_rel_vars)
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
  obs_kurtosis_in_range <- min(sim_kurtosis) < observed_kurtosis && max(sim_kurtosis) > observed_kurtosis

  # We must also account for this when setting plot boundaries and histogram bins
  mean_lower <- ifelse(min(sim_means) < observed_mean, min(sim_means), observed_mean)
  mean_upper <- ifelse(max(sim_means) > observed_mean, max(sim_means), observed_mean)
  if ( obs_mean_in_range ) {
    mean_breaks <- seq(mean_lower,mean_upper, length.out=50)
  } else {
    mean_breaks <- seq(mean_lower,mean_upper, length.out=500)
  }
  
  var_lower <- ifelse(min(sim_vars) < observed_var, min(sim_vars), observed_var)
  var_upper <- ifelse(max(sim_vars) > observed_var, max(sim_vars), observed_var)
  if ( obs_var_in_range ) {
    var_breaks <- seq(var_lower,var_upper, length.out=50)
  } else {
    var_breaks <- seq(var_lower,var_upper, length.out=500)
  }
  
  rel_var_lower <- ifelse(min(sim_rel_vars) < observed_rel_var, min(sim_rel_vars), observed_rel_var)
  rel_var_upper <- ifelse(max(sim_rel_vars) > observed_rel_var, max(sim_rel_vars), observed_rel_var)
  if ( obs_rel_var_in_range ) {
    rel_var_breaks <- seq(rel_var_lower,rel_var_upper, length.out=50)
  } else {
    rel_var_breaks <- seq(rel_var_lower,rel_var_upper, length.out=500)
  }
  
  kurtosis_lower <- ifelse(min(sim_kurtosis) < observed_kurtosis, min(sim_kurtosis), observed_kurtosis)
  kurtosis_upper <- ifelse(max(sim_kurtosis) > observed_kurtosis, max(sim_kurtosis), observed_kurtosis)
  if ( obs_kurtosis_in_range ) {
    kurtosis_breaks <- seq(kurtosis_lower,kurtosis_upper, length.out=25)
  } else {
    kurtosis_breaks <- seq(kurtosis_lower,kurtosis_upper, length.out=250)
  }

  # Arrange a grid of histograms
  par(mfrow=c(2,2))
  
  # Plot the histograms
  hist(sim_means,main="Distribution of test statistic",xlab="mean",breaks=mean_breaks,border=NA,col="deepskyblue2")
  abline(v=observed_mean,col="darkorange2",lwd=1.2)

  hist(sim_vars,main="Distribution of test statistic",xlab="variance",breaks=var_breaks,border=NA,col="deepskyblue2")
  abline(v=observed_var,col="darkorange2",lwd=1.2)
  
  hist(sim_rel_vars,main="Distribution of test statistic",xlab="var/mean",breaks=rel_var_breaks,border=NA,col="deepskyblue2")
  abline(v=observed_rel_var,col="darkorange2",lwd=1.2)
  
  hist(sim_kurtosis,main="Distribution of test statistic",xlab="kurtosis",breaks=kurtosis_breaks,border=NA,col="deepskyblue2")
  abline(v=observed_kurtosis,col="darkorange2",lwd=1.2)
  
  # Undo our messing up of par() in case user wants non-gridded plot later
  par(mfrow=c(1,1))
  
}
