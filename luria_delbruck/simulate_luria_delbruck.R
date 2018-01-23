# This is a function that will simulate according to the Luria-Delbruck model (as well as more general, related models)
# number.of.experiments determines how many total experiments are simulated
# number.of.tubes determines how many tubes were reared in each experiment
# There are a number of sets of parameters to play with
# 1) The termination conditions
#    * tube.volume.microliters: how much liquid is in the test tube the bacteria are grown in?
#    * saturation.concentration.per.mL: the concentration (bacteria per microliter) at which growth stops
#    * simulate.for.fixed.time: the default is to simulate until saturation, this allows
#    * end.time: with simulate.for.fixed.time allows you to stop growth at some number of hours
# 2) Factors affecting the appearence and growth of mutants
#    * sensitive.doubling.time.hours: how many hours to double the size of the wildtype population?
#    * doubling.time.ratio: the doubling time of the resistant divided by the doubling time of the wildtype
#    * mutation.probability: the probability that, when a wiltype divides, it produces a mutant offspring
#    * intial.sensitive.population.size.per.test.tube: how many sensitive bacteria are in the innoculant (assume 0 resistant)?
# 3) Random fluctuations during plating process
#    * enable.tube.to.plate.variability: allows us to model variability induced by pipetting technique
#    * excess.tube.to.plate.standard.deviation: the bigger this value, the more variability (there is variability at 0, more as this term increases)
# 4) Induced mutation (this allows us to include the Poisson model in our simulations)
#    * induced.mutation.probability: the probability that a sensitive bacterium mutates to become resistant, post-plating, because of addition of virus
#    * induced.mutation.standard.deviation: allows the rate of mutation to vary accross individuals, at 0, the probability is constant, the higher the value the more the rate varies around induced.mutation.probability
LuriaDelbruckInSilico <- function(number.of.experiments,
                                  number.of.tubes,
                                  tube.volume.microliters=300,
                                  plated.volume.microliters=100,
                                  saturation.concentration.per.mL=1e9,
                                  mutation.probability=3e-9,
                                  sensitive.doubling.time.hours=1,
                                  doubling.time.ratio=1,
                                  intial.sensitive.population.size.per.test.tube=10,
                                  induced.mutation.probability=0,
                                  induced.mutation.standard.deviation=0,
                                  enable.tube.to.plate.variability=FALSE,
                                  excess.tube.to.plate.standard.deviation=0,
                                  simulate.for.fixed.time=FALSE,
                                  end.time=NA) {
  # recover()
  
  #################
  # Preliminaries #
  #################
  
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
  
  if ( excess.tube.to.plate.standard.deviation != 0 && enable.tube.to.plate.variability == FALSE ) {
    stop("Argument \"excess.tube.to.plate.standard.deviation\" requires argument \"enable.tube.to.plate.variability\" to be TRUE")
  }
  
  if ( simulate.for.fixed.time && (is.na(end.time) || end.time < 0) ) {
    stop("To simulate for fixed time, specify end time > 0")
  }
  
  ###############
  # Useful math #
  ###############
  
  # Find population size limit in test tubes
  maximum_bacteria_per_tube <- (saturation.concentration.per.mL/1000)*tube.volume.microliters

  # Get birth rate from doubling time
  normal_reproduction_rate <- log(2)/sensitive.doubling.time.hours
  
  # When drawing a sample of media to plate, there is variance in how many bacteria are selected
  # Here we model this with a negative binomial distribution
  # This allows us to inflate the variance over the standard Poisson/binomial models
  # We first need the mean (mu) of the distribution, given by the concentration
  # The phi parameter is inversely proportional to this excess variability, so we flip it around
  tube_to_plate_mu <- plated.volume.microliters * saturation.concentration.per.mL/1000
  tube_to_plate_phi <- 1/excess.tube.to.plate.standard.deviation
  
  # Go from ratio of doubling time to ratio of reproduction rates (twice as long to double is haf as fast to reproduce)
  reproduction_rate_ratio <- 1/doubling.time.ratio
  
  ###################################################
  # Simulate a number of Luria-Delbruck experiments #
  ###################################################
  
  # We will store the results of each replicate experiment as a column in a matrix
  n_resistant_final <- matrix(NA,nrow=number.of.tubes,ncol=number.of.experiments)
  
  # recover()
  
  # Each loop is now simulating a single experiment
  for (rep in 1:number.of.experiments) {
    
    #####################
    # Grow the bacteria #
    #####################
    
    starting_bacteria <- matrix(NA,nrow=number.of.tubes,ncol=2)
    
    if ( mutation.probability > 0 ) {  
      # We'll store the simulations in a matrix, first column is number of wild-type, second is mutants
      for (i in 1:number.of.tubes) {
        starting_bacteria[i,] <- simulateBacterialGrowth(final.population.size=maximum_bacteria_per_tube,mutation.probability=mutation.probability,normal.reproduction.rate=normal_reproduction_rate,reproduction.rate.ratio=reproduction_rate_ratio,intial.sensitive.population.size.per.test.tube=intial.sensitive.population.size.per.test.tube,simulate.for.fixed.time=simulate.for.fixed.time,end.time=end.time)
      }
    } else {
      # No mutants, we know how many sensitve bacteria there will be from the growth rate
      if ( simulate.for.fixed.time ) {
        number_at_t_stop <- intial.sensitive.population.size.per.test.tube*exp(end.time*normal_reproduction_rate)
        if ( number_at_t_stop > maximum_bacteria_per_tube) {
          number_sensitive <- maximum_bacteria_per_tube
        } else {
          number_sensitive <- number_at_t_stop
        }
      } else {
        number_sensitive <- maximum_bacteria_per_tube
      }
        
      starting_bacteria[,1] <- number_sensitive
      starting_bacteria[,2] <- 0
    }
    ######################
    # Plate the bacteria #
    ######################
    
    # Draw the number of bacteria/colony founding units that are plated 
    if ( enable.tube.to.plate.variability ) {
      # Under the variable model, we don't get exactly the same fraction every time
      n_cfu <- rnbinom(n=number.of.tubes,mu=tube_to_plate_mu,size=tube_to_plate_phi)
    } else {
      # Here we assume we get exactly plated.volume.microliters/tube.volume.microliters in each pipette
      n_cfu <- tube_to_plate_mu
    }
  
    # The probability that any one of those bacteria is resistant is the frequency of resistant bacteria in that tube
    p_resistant <- starting_bacteria[,2]/rowSums(starting_bacteria)
    
    # Thus the number of resistant bacteria is binomially distributed
    n_resistant_initial <- rbinom(n=number.of.tubes,size=n_cfu,prob=p_resistant)
    
    # Account for directed mutation
    if ( induced.mutation.standard.deviation != 0) {
      # Here our model includes the mutation rate varying from bacterium to bacterium
      # It varies according to a beta(a,b) distribution whose parameters we now compute
      # We exploit the fact that increasing (a+b) decreases the variability in the mutation probability
      # So what we're calling "induced.mutation.standard.deviation" is actually 1/(a+b), and we know the mean we want, and mean = a/(a+b)
      a <- induced.mutation.standard.deviation/induced.mutation.standard.deviation
      b <- (1 - induced.mutation.probability)/induced.mutation.standard.deviation
      mu_induced <- rbeta(number.of.tubes*(n_cfu-n_resistant_initial),a,b)
    } else {
      mu_induced <- induced.mutation.probability
    }
    
    n_resistant_final[,rep] <- n_resistant_initial + rbinom(n=number.of.tubes,size=n_cfu-n_resistant_initial,prob=mu_induced)
    # Does it seem weird that we're not using the Poisson distribution here?
    # The Poisson and the binomial are essentially equivalent when the mutation rate is very small and there are many bacteria
    # But, the Poisson can generate any counting number
    # Since we've already specified a number of susceptible bacteria, it would be weird to allow ourselves to suddenly have more direct-mutation-caused resistant bacteria than we had susceptible
    # A curious person could plug in the Poisson instead and see what happens (the probability of it happening in any one sample is very small, but over many simulated experiments with many simulated tubes, it is higher)
  }
  
  return(n_resistant_final)
}

# This function simulates the growth of bacteria in test tubes
# It is separate from the above function so that the above function is "only" 130 lines, not 250
# Parameters: 
#             final.population.size: the population size at saturation
#             mutation.probability: the probability that a sensitive bateria produces a resistant when dividing
#             normal.reproduction.rate: r in e^rt for exponential growth of bacteria
#             intial.sensitive.population.size.per.test.tube: the starting number of the wild-type bacteria
#             simulate.for.fixed.time: should we stop at a given time, rather than at saturation?
#             end time: the (optional) time at which growth is terminated

simulateBacterialGrowth <- function(final.population.size,
                                    mutation.probability,
                                    normal.reproduction.rate,
                                    reproduction.rate.ratio,
                                    intial.sensitive.population.size.per.test.tube,
                                    simulate.for.fixed.time=FALSE,
                                    end.time=NA) {
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
  t_saturation <- log(final.population.size/intial.sensitive.population.size.per.test.tube)/(B_w)
  
  # Get the times
  mutant_origin_times <- generateLuriaDelbruckPoissonProcess(mu=mutation.probability,N0=intial.sensitive.population.size.per.test.tube,beta=B_w,end.time=t_saturation)
  
  aliveAtT <- function(t_,n_mutations) {
    alive <- round(intial.sensitive.population.size.per.test.tube*exp(B_w*t_)) + sum(round(exp(B_m*(t_ - mutant_origin_times))))
    return(alive)
  }
  
  # Possible early termination 1: are there mutants?
  if ( length(mutant_origin_times) == 0 ) {
    if ( simulate.for.fixed.time ) {
      # Population size at user-defined end-time
      wt_at_t_stop <- round(intial.sensitive.population.size.per.test.tube*exp(B_w*t_))
      wt_at_t_stop <- min(wt_at_t_stop,final.population.size)
    } else {
      # Saturation size
      wt_at_t_stop <- round(intial.sensitive.population.size.per.test.tube*exp(B_w*t_saturation))
    } 
    mutants_at_t_stop <- 0
    bacteria <- c(wt_at_t_stop,mutants_at_t_stop)
    return(bacteria)
  }
  
  # Possible early termination 1: user-defined end time (before saturation)
  if ( simulate.for.fixed.time ) {
    
    # Count the bacteria
    N_t <- aliveAtT(t_)
    
    # We have checked that the user end time does not exceed the maximum saturation time
    # We have to make sure the user-defined end time isn't past the point of saturation
    t_ <- end.time
    
    if ( N_t < final.population.size ) {
      wt_at_t_stop <- round(intial.sensitive.population.size.per.test.tube*exp(B_w*t_))
      mutants_at_t_stop <- sum(round(exp(B_m*(t_ - mutant_origin_times))))
      bacteria <- c(wt_at_t_stop,mutants_at_t_stop)
      return(bacteria)
    } else {
      warning("Specified end point past time of saturation, simulations will be slower.")
    }
  }
  
  
  
  # Now we find the time that this population reaches saturation
  # We do this numerically (we ave the computer search for us), since we can't write a general enough analytical solution
  # To do this, we need to find a time when we were below saturation (a lower bound for searching)
  
  # Index the times of mutation events (these are good places to test for bounds)
  i <- length(mutant_origin_times)
  
  # We first check if the last mutation event could work (good guess if not a user-defined end-time)
  t_ <- mutant_origin_times[i]
  
  # Count the bacteria
  N_t <- aliveAtT(t_)
  
  
  # We may not enter this loop, that's fine (faster, in fact), it means our first guess was good
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
  t_ <- res$minimum
  
  # Find and return number of mutants at stopping time
  # wt_at_t_stop <- round(intial.sensitive.population.size.per.test.tube*exp(B_w*t_))
  # mutants_at_t_stop <- sum(round(exp(B_m*(t_ - mutant_origin_times))))
  mutants_at_t_stop <- round(sum(exp(B_m*(t_ - mutant_origin_times))))
  wt_at_t_stop <- final.population.size - mutants_at_t_stop
  bacteria <- c(wt_at_t_stop,mutants_at_t_stop)
  
  return(bacteria)
}

# This is a helper function for simulating bacterial growth, it simulates the times of mutations
# Generates event times from 0 to end.time from a Luria-Delbruck-esque NHPP
# Parameters: 
#             mu: the mutation rate/probability
#             beta: the growth rate of the wild-type bacteria
#             N0: the starting number of the wild-type bacteria
#             end time: the time at which growth is terminated
generateLuriaDelbruckPoissonProcess <- function(mu,beta,N0,end.time) {
  
  # The integrated rate allows us to generate a number of events that happened in the interval
  integrated_rate <- mu*N0/beta*(exp(beta*end.time) - 1)
  n_events <- rpois(1,integrated_rate)
  
  # Pasupathy's Algorithm #5
  u <- runif(n_events)
  times <- log(u * exp(beta*end.time) - u + 1)/(beta)
  
  return(sort(times))
}

