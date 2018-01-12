trackP <- function(final.population.size,
                                 mutation.probability,
                                 normal.reproduction.rate,
                                 reproduction.rate.ratio,
                                 intial.sensitive.population.size.per.test.tube,
                                 intial.resistant.population.size.per.test.tube) {
  
  # recover()
  p <- numeric(final.population.size-1)
  #########
  # Setup #
  #########
  
  # birth rates in each class, number of wild-types and mutants
  B_w <- normal.reproduction.rate
  B_m <- normal.reproduction.rate * reproduction.rate.ratio
  
  # number of cells currently in each class
  N_t <- c(intial.sensitive.population.size.per.test.tube,intial.resistant.population.size.per.test.tube)
  names(N_t) <- c("number of wild-type","number of mutants")
  
  # track the per-type and total rate of the birth process
  rate_w <- N_t[1]*B_w
  rate_m <- N_t[2]*B_m
  total_rate <- rate_w + rate_m
  
  ############
  # Simulate #
  ############
  
  # Keep colonies growing until saturation population size is achieved
  for (i in 1:(final.population.size-1)) { 
    
    # Calculate probability that event was a birth in wild-type
    p_mutant_increase <- (rate_w * mutation.probability + rate_m) / total_rate
    p[i] <- p_mutant_increase
    # What kind of event happened?
    # Birth in a wild-type or a mutant?
    birth_in_m <- rbinom(n=1,size=1,p=p_mutant_increase)
    if ( birth_in_m ) {
      N_t[2] <- N_t[2] + 1
      rate_m <- rate_m + B_m
      total_rate <- total_rate + B_m
    } else {
      N_t[1] <- N_t[1] + 1
      rate_w <- rate_w + B_w
      total_rate <- total_rate + B_w
    }
    
    # # Time until next event
    # wt <- rexp(1,total_rate)
    # t_elapsed <- t_elapsed + wt
    
  }
  
  # return(N_t)
  return(p)
}

approxGrowBactSilico <- function(final.population.size,
                                 mutation.probability,
                                 normal.reproduction.rate,
                                 reproduction.rate.ratio,
                                 intial.sensitive.population.size.per.test.tube) {
  
  # recover()
  
  #########
  # Setup #
  #########
  
  # birth rates in each class, number of wild-types and mutants
  B_w <- normal.reproduction.rate
  B_m <- normal.reproduction.rate * reproduction.rate.ratio
  
  # number of cells currently in each class
  N_t <- c(intial.sensitive.population.size.per.test.tube,0)
  names(N_t) <- c("number of wild-type","number of mutants")
  
  # track the per-type and total rate of the birth process
  rate_w <- N_t[1]*B_w
  rate_m <- N_t[2]*B_m

  ############
  # Simulate #
  ############
  
  # Number of WT births before first mutant arises
  p_mutant_increase <- (rate_w * mutation.probability + rate_m) / (rate_w + rate_m)
  
  first_interval <- rgeom(1,p_mutant_increase)
  
  
  if ( (first_interval + 1) >= final.population.size ) {
    N_t[1] <- final.population.size
    return(N_t)
  } else {
    
    # Recalculate
    N_t[1] <- N_t[1] + first_interval
    N_t[2] <- N_t[2] + 1
    rate_w <- N_t[1]*B_w
    rate_m <- N_t[2]*B_m
    
    while ( sum(N_t) < final.population.size ) {
      # Control precision of approximation
      thresh <- 1000*max(1,log10(N_t[1])-4)
                    
      p_mutant_increase <- (rate_w * mutation.probability + rate_m) / (rate_w + rate_m)
      next_interval <- rgeom(1,p_mutant_increase)
      if ( next_interval > thresh || is.na(next_interval)) { # It appears that NAs are generated instead of Infinity
        N_t[1] <- N_t[1] + thresh
        rate_w <- N_t[1]*B_w
      } else {
        N_t[1] <- N_t[1] + next_interval
        N_t[2] <- N_t[2] + 1
        rate_w <- N_t[1]*B_w
        rate_m <- N_t[2]*B_m
        
      }
    }
    
    return(N_t)
  }

}
