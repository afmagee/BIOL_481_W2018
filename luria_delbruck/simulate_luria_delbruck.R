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

# This is a function that will simulate according to the Luria-Delbruck model (as well as more general, related models)
# Guide to 
LuriaDelbruckInSilico <- function(number.of.tubes,
                                  tube.volume.microliters=300,
                                  saturation.concentration.microliters=1e9,
                                  mutation.probability=3e-9,
                                  non.resistant.doubling.time.hours=1,
                                  reproduction.rate.ratio=1,
                                  intial.sensitive.population.size.per.test.tube=10,
                                  induced.mutation.probability=0,
                                  induced.mutation.variability=0,
                                  enable.tube.to.plate.variability=FALSE,
                                  excess.tube.to.plate.variability=0,
                                  verbose=TRUE,
                                  ...) {
  recover()
  
  #################
  # Preliminaries #
  #################
  
  if ( hasArg(intial.resistant.population.size.per.test.tube) ) {
    intial_resistant_population_size_per_test_tube <- list(...)$intial.resistant.population.size.per.test.tube
  } else {
    intial_resistant_population_size_per_test_tube <- 0
  }
  
  if ( hasArg(plated.volume.microliters) ) {
    plated_volume_microliters <- list(...)$plated.volume.microliters
  } else {
    plated_volume_microliters <- 100
  }
  
  # Error checking
  if ( induced.mutation.probability > 1) {
    stop("Probabilities are in [0,1], invalid input to argument \"induced.mutation.probability\" ")
  }
  
  if ( mutation.probability > 1) {
    stop("Probabilities are in [0,1], invalid input to argument \"mutation.probability\" ")
  }
  
  if ( !(enable.tube.to.plate.variability == TRUE || enable.tube.to.plate.variability == FALSE) ) {
    stop("Argument \"enable.tube.to.plate.variability\" must be TRUE or FALSE")
  }
  
  if ( excess.tube.to.plate.variability != 0 && enable.tube.to.plate.variability == FALSE ) {
    stop("Argument \"excess.tube.to.plate.variability\" requires argument \"enable.tube.to.plate.variability\" to be TRUE")
  }
  
  ###############
  # Useful math #
  ###############
  
  # Find population size limit in test tubes
  maximum_bacteria_per_tube <- saturation.concentration.microliters*tube.volume.microliters
  
  # Get birth rate from doubling time
  normal_reproduction_rate <- log(2)/non.resistant.doubling.time.hours
  
  # When drawing a sample of media to plate, there is variance in how many bacteria are selected
  # Here we model this with a negative binomial distribution
  # This allows us to inflate the variance over the standard Poisson/binomial models
  # We first need the mean (mu) of the distribution, given by the concentration
  # The phi parameter is inversely proportional to this excess variability, so we flip it around
  tube_to_plate_mu <- plated_volume_microliters * saturation.concentration.microliters
  tube_to_plate_phi <- 1/excess.tube.to.plate.variability

  #####################
  # Grow the bacteria #
  #####################
  
  starting_bacteria <- matrix(NA,nrow=number.of.tubes,ncol=2)
  for (i in 1:number.of.tubes) {
    starting_bacteria[i,2] <- growBacteriaInSilico(final.population.size=maximum_bacteria_per_tube,mutation.probability=mutation.probability,normal.reproduction.rate=normal.reproduction.rate,reproduction.rate.ratio=reproduction.rate.ratio,intial.sensitive.population.size.per.test.tube=intial.sensitive.population.size.per.test.tube,intial.resistant.population.size.per.test.tube=)
  }

  ######################
  # Plate the bacteria #
  ######################
  
  # Draw number of bacteria/colony founding units that are plated 
  if ( excess_tube_to_plate_variability ) {
    n_cfu <- tube_to_plate_mu
  } else {
    n_cfu <- rnbinom(n=number.of.tubes,mu=tube_to_plate_mu,size=tube_to_plate_phi)
  }

  # The probability that any one of those bacteria is resistant is the frequency of resistant bacteria in that tube
  p_resistant <- starting_bacteria[,2]/maximum_bacteria_per_tube
  
  # Thus the number of resistant bacteria is binomially distributed
  n_resistant_initial <- rbinom(n=number.of.tubes,size=n_cfu,prob=p_resistant)
  
  # Account for directed mutation
  if ( induced.mutation.variability != 0) {
    # Here our model includes the mutation rate varying from bacterium to bacterium
    # It varies according to a beta(a,b) distribution whose parameters we now compute
    # We exploit the fact that increasing (a+b) decreases the variability in the mutation probability
    # So what we're calling "induced.mutation.variability" is actually 1/(a+b), and we know the mean we want, and mean = a/(a+b)
    a <- induced.mutation.variability/induced.mutation.variability
    b <- (1 - induced.mutation.probability)/induced.mutation.variability
    mu_induced <- rbeta(number.of.tubes*(n_cfu-n_resistant_initial),a,b)
  } else {
    mu_induced <- induced.mutation.probability
  }
  
  n_resistant_final <- n_resistant_initial + rbinom(n=number.of.tubes,size=n_cfu-n_resistant_initial,prob=mu_induced)
  # Does it seem weird that we're not using the Poisson distribution here?
  # The Poisson and the binomial are essentially equivalent when the mutation rate is very small and there are many bacteria
  # But, the Poisson can generate any counting number
  # Since we've already specified a number of susceptible bacteria, it would be weird to allow ourselves to suddenly have more direct-mutation-caused resistant bacteria than we had susceptible
  # A curious person could plug in the Poisson instead and see what happens (the probability of it happening in any one sample is very small, but over many simulated experiments with many simulated tubes, it is higher)
  
  return(n_resistant_final)
}

generateLuriaDelbruckPoissonProcess <- function(mu,beta,end.time) {

  # The integrated rate allows us to generate a number of events that happened in the interval
  integrated_rate <- mu/beta*(exp(beta*end.time) - 1)
  n_events <- rpois(1,integrated_rate)
  
  # Pasupathy's Algorithm #5
  u <- runif(n_events)
  
  times <- log(u)/beta + end.time
  
  return(sort(times))
}

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
  
  aliveAtT <- function(t_,n_mutations) {
    alive <- round(intial.sensitive.population.size.per.test.tube*exp(B_w*t_)) + sum(round(exp(B_m*(t_ - mutant_origin_times))))
    return(alive)
  }
  
  # We know at the t_saturation we have too many bacteria (we designed it so)
  # Work backwards until we find a time where we don't have too many bacteria
  i <- length(mutant_origin_times)
  t_ <- mutant_origin_times[i]
  
  # Count the bacteria
  N_t <- aliveAtT(t_)
  
  # We may not enter this loop, that's fine
  while ( N_t > final.population.size ) {
    i <- i - 1
    t_ <- mutant_origin_times[i]
    
    # How many bacteria are alive now?
    N_t <- aliveAtT(t_)

  }
  
  # Now we find the time that this population reaches saturation
  # We do this numerically, since we can't write a general enough analytical solution
  
  # This function is minimized when t_stop is correct
  fn <- function(t_stop) {
    N_t <- aliveAtT(t_stop)
    return(((log(N_t)-log(final.population.size))^2))
  }
  
  # Find the stopping time
  res <- optimize(fn,c(t_,t_saturation))
  
  # Find and return number of mutants at stopping time
  mutants_at_t_stop <- sum(round(exp(B_m*(t_ - mutant_origin_times))))
  
  return(mutants_at_t_stop)
}


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
  # To ensure that we have enough event times, we generate times between the start of the process and the time the wild-types reach saturation
  
  
  # Get the times
  mutant_origin_times <- generateLuriaDelbruckPoissonProcess(mutation.probability,B_w,end.time)
  
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
