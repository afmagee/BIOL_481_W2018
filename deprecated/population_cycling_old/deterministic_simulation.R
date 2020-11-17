H0 <- 10
P0 <- 2

nsteps <- 100000

# These are the values in Bartlet (1960, pg 35)
alpha1 <- 1
beta1 <- 0.1
alpha2 <- 0.5
beta2 <- 0.02

Ht <- numeric(nsteps)
Pt <- numeric(nsteps)

Ht[1] <- H0
Pt[1] <- P0

# If this gets too large the series will start to tend towards infinity
dt_ <- 1e-3

for (t in 1:(nsteps-1)) {
  Ht[t+1] <- Ht[t]+(alpha1 - beta1*Pt[t])*Ht[t]*dt_
  Pt[t+1] <- Pt[t]+(-alpha2 + beta2*Ht[t])*Pt[t]*dt_
}

plot(Ht/Ht[1],type="l")
lines(Pt/Pt[1],lty=3)