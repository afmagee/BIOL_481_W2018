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
#    * resistant.doubling.time.hours: how many hours to double the size of the mutant population?
#    * mutation.probability: the probability that, when a wiltype divides, it produces a mutant offspring
#    * intial.sensitive.population.size.per.test.tube: how many sensitive bacteria are in the innoculant (assume 0 resistant)?
# 3) Random fluctuations during plating process
#    * enable.tube.to.plate.variability: allows us to model variability induced by pipetting technique
# 4) Induced mutation (this allows us to include the Poisson model in our simulations)
#    * induced.mutation.probability: the probability that a sensitive bacterium mutates to become resistant, post-plating, because of addition of virus
#    * induced.mutation.probability.standard.deviation: allows the rate of mutation to vary accross individuals in the population from which the innoculant was taken, thus giving some plates higher mutation rates than others
LuriaDelbruckInSilico <- function(number.of.experiments,
                                  number.of.tubes,
                                  tube.volume.microliters=200,
                                  plated.volume.microliters=30,
                                  saturation.concentration.per.mL=1e9,
                                  mutation.probability=3e-10,
                                  sensitive.doubling.time.hours=1,
                                  resistant.doubling.time.hours=2,
                                  intial.sensitive.population.size.per.test.tube=1,
                                  induced.mutation.probability=1e-8,
                                  induced.mutation.probability.standard.deviation=0,
                                  enable.tube.to.plate.variability=FALSE,
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
  
  # When drawing a sample of media to plate, there may be variance in how many bacteria are selected
  # We would expect to get the number of bacteria in the colony (given by saturation constant and volume)
  tube_to_plate_mean <- plated.volume.microliters * saturation.concentration.per.mL/1000

  # Go from ratio of doubling time to ratio of reproduction rates (twice as long to double is haf as fast to reproduce)
  
  # turn ratio into just "resistant doubling time" 
  mutant_reproduction_rate <- log(2)/resistant.doubling.time.hours
  
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
        starting_bacteria[i,] <- simulateBacterialGrowth(final.population.size=maximum_bacteria_per_tube,mutation.probability=mutation.probability,normal.reproduction.rate=normal_reproduction_rate,mutant.reproduction.rate=mutant_reproduction_rate,intial.sensitive.population.size.per.test.tube=intial.sensitive.population.size.per.test.tube,simulate.for.fixed.time=simulate.for.fixed.time,end.time=end.time)
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
      # It may seem a little weird to allow there to be more than the "mean" (which is the concentration * the volume)
      # One could see this as accounting for minor fluctuations around the saturation constant (or as a modeling simplification)
      n_cfu <- rpois(n=number.of.tubes,lambda=tube_to_plate_mean)
    } else {
      # Here we assume we get exactly plated.volume.microliters/tube.volume.microliters in each pipette
      n_cfu <- tube_to_plate_mean
    }
  
    # The probability that any one of those bacteria is resistant is the frequency of resistant bacteria in that tube
    p_resistant <- starting_bacteria[,2]/rowSums(starting_bacteria)
    
    # Thus the number of resistant bacteria is binomially distributed
    n_resistant_initial <- rbinom(n=number.of.tubes,size=n_cfu,prob=p_resistant)
    
    # Account for directed mutation
    if ( induced.mutation.probability.standard.deviation != 0 ) {
      # Here we have more variation than a simple Poisson distribution
      # We do that by simulating from the Negative Binomial distribution (in two steps, as can be seen as a Poisson with a variable rate, where the rate varies by a gamma distribution)
      # First we get the parameters of the gamma distribution (from mean = mutation probability, and user-defined standard deviation)
      induced_mutation_scale <- (induced.mutation.probability.standard.deviation^2)/induced.mutation.probability
      induced_mutation_shape <- induced.mutation.probability/induced_mutation_scale
      # Now we generate a range of lambdas
      induced_mutation_lambda <- rgamma(number.of.tubes,shape=induced_mutation_shape,scale=induced_mutation_scale) * (n_cfu-n_resistant_initial)
    } else {
      induced_mutation_lambda <- induced.mutation.probability * (n_cfu-n_resistant_initial)
    }
    n_resistant_final[,rep] <- n_resistant_initial + rpois(n=number.of.tubes,lambda=induced_mutation_lambda)
  }
  
  # recover()
  
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
                                    mutant.reproduction.rate,
                                    intial.sensitive.population.size.per.test.tube,
                                    simulate.for.fixed.time=FALSE,
                                    end.time=NA) {
  # recover()
  
  #########
  # Setup #
  #########
  
  # birth rates in each class, number of wild-types and mutants
  B_w <- normal.reproduction.rate
  B_m <- mutant.reproduction.rate
  
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

# For those who explored the directed mutation plus variable mutation rate model, this function may help
# It shows you how likely a range of mutation rates are under the settings of this model
# The plotted dashed green line is the average
# The inputs induced.mutation.probability and induced.mutation.probability.standard.deviation work as in LuriaDelbruckInSilico
# evaluate.probabilities.at: function will report the probability that a bacterium has a mutation rate greater than these values
visualizeInducedMutationVariance <- function(induced.mutation.probability,induced.mutation.probability.standard.deviation,evaluate.probabilities.at=c(1e-10,1e-9,1e-8,1e-7,1e-6)) {
  # recover()
  
  # Gamma parameters
  induced_mutation_scale <- (induced.mutation.probability.standard.deviation^2)/induced.mutation.probability
  induced_mutation_shape <- induced.mutation.probability/induced_mutation_scale
  
  # Find an upper limit to the plot
  dens <- dgamma(seq(0,1e-4,1e-4/1000),shape=induced_mutation_shape,scale=induced_mutation_scale)
  if ( dens[1000] > 1e-4 ) {
    upper <- 1e-4
  } else {
    upper <- min(which(dens < 1e-4))*1e-4/1000
  }
  
  cat("Probability of a bacterium having a mutation rate above specified values:\n")
  for (i in 1:length(evaluate.probabilities.at)) {
    cat(evaluate.probabilities.at[i],": ",pgamma(evaluate.probabilities.at[i],shape=induced_mutation_shape,scale=induced_mutation_scale,lower.tail=FALSE),"\n",sep="")
  }
      
  
  ggplot(data.frame(x=seq(0,upper,length.out=10001),y=dgamma(seq(0,upper,length.out=10001),shape=induced_mutation_shape,scale=induced_mutation_scale))) + 
    geom_line(aes(x,y)) + 
    ggtitle("Distribution of mutation rates among bacteria") + 
    xlab("mutation rate") + 
    ylab("probability density") + 
    geom_vline(xintercept=induced.mutation.probability,color="springgreen3",lwd=1.0,lty=2)
}

