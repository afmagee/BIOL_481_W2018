# This is a function that will simulate according to the Luria-Delbruck model (as well as more general, related models)

# birth_prey <- 0.1
# birth_pred <- 0.1
# 
# death_prey <- 0.1
# death_pred <- 0.1

beta_ratio <- 0.1
beta_1 <- 0.1

pops <- matrix(nrow=1,ncol=3)
pops[1,] <- c(25,2,0)

current <- pops[1,]

simT <- 1000

t_ <- 0

while (t_ < simT) {
  
  
  H_up_P_static   <- current[1]
  H_down_P_static <- (beta_1-beta_1*beta_ratio)*prod(current[1:2]) 
  H_down_P_up     <- beta_1*beta_ratio*prod(current[1:2])
  H_static_P_down <- 0.5*current[2]
  
  wt <- rexp(1,H_up_P_static+H_down_P_static+H_down_P_up+H_static_P_down)
  
  type <- sample(1:4,1,prob=c(H_up_P_static,H_down_P_static,H_down_P_up,H_static_P_down))
  
  if ( type == 1) {
    current[1] <- current[1] + 1
  } else if ( type == 2) {
    current[1] <- current[1] - 1
  } else if ( type == 3) {
    current[2] <- current[2] + 1
  } else if ( type == 4) {
    current[2] <- current[2] - 1
  }
  
  t_ <- t_ + wt
  
  current[3] <- t_
  
  pops <- rbind(pops,current)
  
  if ( prod(current[1:2]) == 0 ) {
    break
  }
  
}



plot(pops[,3],pops[,1],xlab="time",ylab="pop",type="l",ylim=c(0,max(pops[,-3])))
lines(pops[,3],pops[,2],xlab="time",ylab="pop",type="l",lty=2)


plot(pops[,1],pops[,2],type="l")


predPreyTimeseriesLnL <- function(pred,prey,times,TPM) {
  
}


