simulatePopulationCycling <- function(H0=10,
                                      P0=2,
                                      alpha1=1,
                                      beta1=0.1,
                                      alpha2=0.5,
                                      beta2=0.02,
                                      noise.level=2,
                                      elapsed.time=10,
                                      num.replicates=100) {
  
  # recover()
  
  # If this gets too large the series will start to tend towards infinity
  dt_ <- 1e-3
  
  # Storing all the datapoints is onerous and plotting them is intensive
  thin_to <- 1000
  
  nsteps <- elapsed.time/dt_
  
  # Variance of the stohastic part of the stochastic DE
  sig <- noise.level*dt_
  
  # Get the realizations
  prey <- as.data.frame(matrix(nrow=thin_to,ncol=num.replicates))
  pred <- as.data.frame(matrix(nrow=thin_to,ncol=num.replicates))
  
  pb <- txtProgressBar(min = 0, max = num.replicates, style = 3)
  for (i in 1:num.replicates) {
    Ht <- numeric(nsteps)
    Pt <- numeric(nsteps)
    
    Ht[1] <- H0
    Pt[1] <- P0
    
    # Save time by batch-computing noise
    hnoise <- sig*rnorm(nsteps-1,0,sig)
    pnoise <- sig*rnorm(nsteps-1,0,sig)
    
    # Simulate
    for (t in 1:(nsteps-1)) {
      Ht[t+1] <- Ht[t]+(alpha1 - beta1*Pt[t])*Ht[t]*dt_ + hnoise[t]
      Pt[t+1] <- Pt[t]+(-alpha2 + beta2*Ht[t])*Pt[t]*dt_ + pnoise[t]
    }
    
    # Thin out the realizations of the stochastic process
    prey[,i] <- Ht[seq(1,nsteps,length.out=thin_to)]
    pred[,i] <- Pt[seq(1,nsteps,length.out=thin_to)]
    
    setTxtProgressBar(pb, i)
    
  }
  # recover()
  # Get the expectation
  Ht <- numeric(nsteps)
  Pt <- numeric(nsteps)
  
  Ht[1] <- H0
  Pt[1] <- P0
  
  
  for (t in 1:(nsteps-1)) {
    Ht[t+1] <- Ht[t]+(alpha1 - beta1*Pt[t])*Ht[t]*dt_
    Pt[t+1] <- Pt[t]+(-alpha2 + beta2*Ht[t])*Pt[t]*dt_
  }
  
  Ht <- Ht[seq(1,nsteps,length.out=thin_to)]
  Pt <- Pt[seq(1,nsteps,length.out=thin_to)]
  
  plot(x=seq(1,elapsed.time,length.out=thin_to),y=Ht,type="l",col="blue",ylab="population size(s)",xlab="time",ylim=c(0,max(prey)))
  lines(x=seq(1,elapsed.time,length.out=thin_to),y=Pt,type="l",col="red")
  
  matplot(x=seq(1,elapsed.time,length.out=thin_to),y=prey,type="l",col="#0000ff50",lty=1,lwd=0.5,add=TRUE)
  matplot(x=seq(1,elapsed.time,length.out=thin_to),y=pred,type="l",col="#ff000050",lty=1,lwd=0.5,add=TRUE)
  
  # plot(Ht/Ht[1],type="l")
  # lines(Pt/Pt[1],lty=3)
  # abline(h=0)
}

system.time(simulatePopulationCycling(elapsed.time=50,num.replicates=20))



simulatePopulationCycling(elapsed.time=50,num.replicates=20,noise.level=20)


simulatePopulationCycling(elapsed.time=50,num.replicates=100,noise.level=20)


lynxhare <- read_csv("~/Downloads/lynxhare.csv")
lynxhare <- lynxhare[,1:3]
lynxhare <- dplyr::filter(lynxhare, complete.cases(lynxhare))
lynxhare <- dplyr::filter(lynxhare, year > 1896)

plot(lynxhare$year,lynxhare$hare/mean(lynxhare$hare),col="blue",pch=16)
points(lynxhare$year,lynxhare$lynx/mean(lynxhare$lynx),col="red",pch=16)


fitPopulationCycle <- function(time.series) {
  # If this gets too large the series will start to tend towards infinity
  dt_ <- 1e-3
  
  # Storing all the datapoints is onerous and plotting them is intensive
  thin_to <- dim(time.series)[1]
  
  nsteps <- dim(time.series)[1]/dt_
  
  E_Ht <- numeric(nsteps)
  E_Pt <- numeric(nsteps)
  
  E_Ht[1] <- H0
  E_Pt[1] <- P0
  
  
  for (t in 1:(nsteps-1)) {
    Ht[t+1] <- Ht[t]+(alpha1 - beta1*Pt[t])*Ht[t]*dt_
    Pt[t+1] <- Pt[t]+(-alpha2 + beta2*Ht[t])*Pt[t]*dt_
  }
  
  Ht <- Ht[seq(1,nsteps,length.out=thin_to)]
  Pt <- Pt[seq(1,nsteps,length.out=thin_to)]
  
  plot(x=seq(1,elapsed.time,length.out=thin_to),y=Ht,type="l",col="blue",ylab="population size(s)",xlab="time",ylim=c(0,max(prey)))
  lines(x=seq(1,elapsed.time,length.out=thin_to),y=Pt,type="l",col="red")
  
  matplot(x=seq(1,elapsed.time,length.out=thin_to),y=prey,type="l",col="#0000ff50",lty=1,lwd=0.5,add=TRUE)
  matplot(x=seq(1,elapsed.time,length.out=thin_to),y=pred,type="l",col="#ff000050",lty=1,lwd=0.5,add=TRUE)
  
}