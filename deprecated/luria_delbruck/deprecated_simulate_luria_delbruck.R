simulateFixedTimeLD <- function(end.time,
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
  
  # Now, we need the times at which mutations occur
  # These times come from a Non-homogenous Poisson process with a rate equal to the mutation rate times the growth rate of the wild-types
  
  # Get the times
  mutant_origin_times <- generateLuriaDelbruckPoissonProcess(mutation.probability,B_w,intial.sensitive.population.size.per.test.tube,end.time)
  
  mutants_at_t_stop <- sum(exp(B_m*(end.time - mutant_origin_times)))
  
  return(mutants_at_t_stop)
}


calculateLDmoments <- function(end.time,
                               mutation.probability,
                               normal.reproduction.rate,
                               reproduction.rate.ratio,
                               intial.sensitive.population.size.per.test.tube) {
  
  B_w <- normal.reproduction.rate
  B_m <- normal.reproduction.rate * reproduction.rate.ratio
  
  if ( reproduction.rate.ratio == 1 ) {
    E_Xt <- mutation.probability*end.time*exp(B_w*end.time)
    Var_Xt <- mutation.probability/B_w * exp(B_w*end.time) * (exp(B_w*end.time) - 1)
  }
  
  res <- c(E_Xt,Var_Xt)
  names(res) <- c("E[X(t)]","Var[X(t)]")
  
  return(res)
  
}

# Function to simulate bacterial growth, takes the following arguments
#   final.population.size: the population size (non-resistant + resistant) at which we stop simulating
#   mutation.probability: the probability that, in dividing, a non-resistant produces a resistant offspring
#   normal.reproduction.rate: reproduction rate of the non-resistant, as a rate for a birth-death process
#   reproduction.rate.ratio: rate(non-resistant)/rate(resistant), ratio means we can think of it in fitness terms
#   intial.sensitive.population.size.per.test.tube: as named
#   intial.resistant.population.size.per.test.tube: as named

growBacteriaInSilico <- function(final.population.size,
                                 mutation.probability,
                                 normal.reproduction.rate,
                                 reproduction.rate.ratio,
                                 intial.sensitive.population.size.per.test.tube,
                                 intial.resistant.population.size.per.test.tube) {
  
  # recover()
  
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
  
  return(N_t)
  
}

# To simulate the whole experiment, we must be able to simulate the growth of the bacteria, which we do here

# We have B_m the mutant reproduction time, B_w the wild-type reproduction time
# To recover the classical Luria-Delbruck, we assume that r = B_m/B_w = 1
# We also have N0 the number of starting cells
# We get original L-D by assuming N0=1
# Lastly we have p, the probability that a wild-type produces one offspring of each type

simulateBacteriaSlow <- function(total.time,
                                 mutation.probability,
                                 normal.reproduction.rate,
                                 reproduction.rate.ratio=1,
                                 intial.population.size.per.test.tube=1,
                                 return.only.last.pop.sizes=TRUE) {
  
  # recover()
  
  #########
  # Setup #
  #########
  
  # birth rates in each class, number of wild-types and mutants
  B_w <- normal.reproduction.rate
  B_m <- normal.reproduction.rate * reproduction.rate.ratio
  
  # number of cells currently in each class
  N_w <- intial.population.size.per.test.tube
  N_m <- 0
  
  # track the total rate of the birth process
  total_rate <- N_w*B_w + N_m*B_m
  
  # Track the populations over time, starting at t0
  N_t <- matrix(c(N_w,N_m),nrow=1,ncol=2)
  colnames(N_t) <- c("number of wild-type","number of mutants")
  
  # Track the event times
  t_elapsed <- 0
  event_times <- 0
  
  ############
  # Simulate #
  ############
  
  # Draw a waiting time
  wt <- rexp(1,total_rate)
  t_elapsed <- t_elapsed + wt
  
  
  # If the next event would carry us past the final time of the experiment, we're done
  # If this were any other language, this would be a do-while loop, without a simulated time before the loop even began
  while (t_elapsed < total.time) { 
    
    # Track time
    event_times <- c(event_times,t_elapsed)
    
    # Calculate probability that event was a birth in wild-type
    p_wildtype <- (B_w * N_w) / (B_w * N_w + B_m * N_m)

    # What kind of event happened?
    birth_in_w <- rbinom(n=1,size=1,p=p_wildtype)
    if ( birth_in_w ) {
      was_mutation <- rbinom(n=1,size=1,p=mutation.probability)
      if ( was_mutation ) {
        N_m <- N_m + 1
        total_rate <- total_rate + B_m
      } else {
        N_w <- N_w + 1
        total_rate <- total_rate + B_w
      }
    } else {
      N_m <- N_m + 1
      total_rate <- total_rate + B_m
    }
    
    N_t <- rbind(N_t,c(N_w,N_m))
    
    # Time until next event
    wt <- rexp(1,total_rate)
    t_elapsed <- t_elapsed + wt
    
  }
  
  if ( return.only.last.pop.sizes ) {
    N <- c(N_w,N_m)
    names(N) <- c("number of wild-type","number of mutants")
    return(N)
  } else {
    return(list(N.of.t=N_t,event.times=event_times))
  }
  
}

simulateBacteriaFast <- function(total.time,
                                 mutation.probability,
                                 normal.reproduction.rate,
                                 reproduction.rate.ratio=1,
                                 intial.population.size.per.test.tube=1) {
  # The slow way is troubled by a few issues
  # 1) storage: exponential growth requires storing a ton of times and population sizes
  # 2) inefficiency of constant concatenating
  #  The solution: kill 2 birds with one stone, store only state at current time and current time
  # 3) repeated calls to rexp are wasteful
  #  Solution less obvious, as any way around this requires additional loops/loop statements
  # 4) repeated math is wasteful
  #  Solution: minimize waste
  
  # Current implementation speed gains: 10-fold (or greater)
  
  # recover()
  
  #########
  # Setup #
  #########
  
  # birth rates in each class, number of wild-types and mutants
  B_w <- normal.reproduction.rate
  B_m <- normal.reproduction.rate * reproduction.rate.ratio
  
  # number of cells currently in each class
  N_t <- c(intial.population.size.per.test.tube,0)
  names(N_t) <- c("number of wild-type","number of mutants")
  
  # track the per-type and total rate of the birth process
  rate_w <- N_t[1]*B_w
  rate_m <- N_t[2]*B_m
  total_rate <- rate_w + rate_m
  
  # Track the event times
  t_elapsed <- 0

  ############
  # Simulate #
  ############
  
  # Draw a waiting time
  wt <- rexp(1,total_rate)
  t_elapsed <- t_elapsed + wt
  
  # If the next event would carry us past the final time of the experiment, we're done
  # If this were any other language, this would be a do-while loop, without a simulated time before the loop even began
  while (t_elapsed < total.time) { 
    
    # Calculate probability that event was a birth in wild-type
    p_wildtype <- rate_w / total_rate
    
    # What kind of event happened?
    birth_in_w <- rbinom(n=1,size=1,p=p_wildtype) # still faster than a uniform+ifelse
    if ( birth_in_w ) {
      was_mutation <- rbinom(n=1,size=1,p=mutation.probability)
      if ( was_mutation ) {
        N_t[2] <- N_t[2] + 1
        rate_m <- rate_m + B_m
        total_rate <- total_rate + B_m
      } else {
        N_t[1] <- N_t[1] + 1
        rate_w <- rate_w + B_w
        total_rate <- total_rate + B_w
      }
    } else {
      N_t[2] <- N_t[2] + 1
      rate_m <- rate_m + B_m
      total_rate <- total_rate + B_m
    }
    
    # Time until next event
    wt <- rexp(1,total_rate)
    t_elapsed <- t_elapsed + wt
    
  }
  
  return(N_t)
  
}

# This is a function that will simulate according to the Luria-Delbruck model (as well as more general, related models)
simulateSingleTubeLuriaDelbruck <- function(number.of.plates,
                                            total.growth.time,
                                            mutation.probability,
                                            normal.reproduction.rate,
                                            reproduction.rate.ratio=1,
                                            number.of.colonies.per.plate=function(){return(20)},
                                            number.of.founding.bacteria.per.colony=function(){return(1)},
                                            intial.population.size.per.test.tube=1) {
  # recover()
  
  # First we need to grow the bacteria in the test-tube
  starting_bacteria <- simulateBacteriaFast(total.time=total.growth.time,mutation.probability=mutation.probability,normal.reproduction.rate=normal.reproduction.rate,reproduction.rate.ratio=reproduction.rate.ratio,intial.population.size.per.test.tube=intial.population.size.per.test.tube)
  cat(starting_bacteria,"\n")
  # Now we plate them and expose them to the virus, then record number of surviving colonies
  n_surviving_colonies <- numeric(number.of.plates)
  for (i in 1:number.of.plates) {
    number_of_colonies_on_this_plate <- number.of.colonies.per.plate()
    number_of_mutants <- numeric(number_of_colonies_on_this_plate)
    # What is the probability of getting a mutant? It's just the fraction of mutants
    p_mutant <- starting_bacteria[2]/(starting_bacteria[1]+starting_bacteria[2])
    for (j in 1:number_of_colonies_on_this_plate) {
      # Sample the appropriate number of bacteria
      how_many_founders <- number.of.founding.bacteria.per.colony()
      number_of_mutants[j] <- rbinom(n=1,size=how_many_founders,p=p_mutant)
    }
    n_surviving_colonies[i] <- sum(number_of_mutants > 0)
  }
  return(n_surviving_colonies)
}

s <- numeric(10)
for (i in 1:10) {
  s[i] <- simulateSingleTubeLuriaDelbruck(1,5,0.05,1,number.of.colonies.per.plate=function(){return(rpois(1,100))},number.of.founding.bacteria.per.colony=function(){return(rpois(1,5))})
}

var(s)
mean(s)

simulateBacteriaSlow(total.time=10,mutation.probability=0.1,normal.reproduction.rate=1,reproduction.rate.ratio=1)

simulateBacteriaFast(total.time=10,mutation.probability=0.1,normal.reproduction.rate=1,reproduction.rate.ratio=1)


system.time({
  for (i in 1:20) {
    simulateBacteriaSlow(total.time=10,mutation.probability=0.1,normal.reproduction.rate=1,reproduction.rate.ratio=1)
  }
})

system.time({
  for (i in 1:20) {
    simulateBacteriaFast(total.time=10,mutation.probability=0.1,normal.reproduction.rate=1,reproduction.rate.ratio=1)
  }
})


# This is probably better cast as a birth/birth-death process (where we fix the death rate to be 0) for simulating in continuous time)
# Better at least in that we could perform inference under this model

simulateDeterministicStochastic <- function(final.population.size,
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
  
  # Now, we need the times at which mutations occur
  # These times come from a Non-homogenous Poisson process with a rate equal to the mutation rate times the growth rate of the wild-types
  # To ensure that we have enough event times, we generate times between the start of the process and the time the wild-types reach saturation
  
  
  #### FIX THIS this assumes initially 1 wildtype (and 0 mutants)
  t_saturation <- log(final.population.size)/B_w
  
  # Get the times
  mutant_origin_times <- generateLuriaDelbruckPoissonProcess(mutation.probability,B_w,t_saturation)
  
  # Make a list of functions that will tell us how many bacteria are in each lineage
  lineages <- list()
  
  # The main lineage, non-resistant bacteria
  lineages[[1]] <- function(t_) {intial.sensitive.population.size.per.test.tube*exp(B_w*t_)}
  
  # The mutants, which grow starting at t - time_that_mutant_arose
  for (i in 1+(1:length(mutant_origin_times))) {
    lineages[[i]] <- function(t_,index){exp(B_m*(t_ - mutant_origin_times[index]))}
  }
  
  # We know at the t_saturation we have too many bacteria (we designed it so)
  # Work backwards until we find a time where we don't have too many bacteria
  i <- length(mutant_origin_times)
  t_ <- mutant_origin_times[i]
  
  # Count the bacteria
  N_t <- lineages[[1]](t_)
  for (j in 2:i) {
    N_t <- N_t + lineages[[j]](t_,j-1)
  }
  
  while ( N_t > final.population.size ) {
    i <- i - 1
    t_ <- mutant_origin_times[i]
    
    # How many bacteria are alive now?
    N_t <- lineages[[1]](t_)
    for (j in 2:i) {
      N_t <- N_t + lineages[[j]](t_,j-1)
    }
    
  }
  
  # Now we find the time that this population reaches saturation
  # We do this numerically, since we can't write a general enough analytical solution
  
  # Make a function to calculate how far off we'd be stopping at time t_stop
  
  # 
  #   # This function is 0 when t_stop is the time at which there are final.population.size bacteria 
  #   fn <- function(t_stop) {
  #     N_t <- lineages[[1]](t_stop)
  #     for (j in 2:i) {
  #       N_t <- N_t + lineages[[j]](t_stop,j-1)
  #     }
  #     return(((N_t/final.population.size)^2-1))
  #   }
  #   
  #   # Find the stopping time
  #   val <- uniroot(fn,c(t_,t_saturation))
  
  # This function is minimized when t_stop is correct
  fn <- function(t_stop) {
    N_t <- lineages[[1]](t_stop)
    for (j in 2:i) {
      N_t <- N_t + lineages[[j]](t_stop,j-1)
    }
    return(((log(N_t)-log(final.population.size))^2))
  }
  
  # Find the stopping time
  res <- optimize(fn,c(t_,t_saturation))
  
  # Find and return number of mutants at stopping time
  N_t <- 0
  for (j in 2:i) {
    N_t <- N_t + lineages[[j]](res$minimum,j-1)
  }
  
  return(N_t)
}

# too slow
generateLuriaDelbruckPoissonProcess <- function(mu,beta,end.time) {
  # There's no need to compute these numbers thousands of times
  beta_over_mu <- beta/mu
  beta_inverse <- 1/beta

  # Pasupathy's Algorithm #3 (Cinlar’s Method)
  i <- 0
  s  <- 0
  t_ <- 0
  times <- c()
  while ( t_ < end.time ) {
    i <- i + 1
    s <- s - log(runif(1))
    # This is the thing we're computing, written clearly
    # times[i] <- 1/beta * log(s * beta/mu + 1)
    # This is the thing we're computing, written efficiently
    t_ <-  beta_inverse * log(s * beta_over_mu + 1)
    times[i] <- t_
  }
  return(times)
}