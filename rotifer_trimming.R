ar.orig <- read.csv("~/Downloads/Default Dataset (1).csv",header=FALSE)
plot(ar.orig)

al <- ar.orig[1:39,]
ro <- ar.orig[40:78,]

cbind(al[,1],ro[,1])

setdiff(round(al[,1]),25:63)

ar.data <- data.frame(times=25:63,prey=al[,2],pred=ro[,2])

ar.data[,2:3] <- ar.data[,2:3]*380

matplot(ar.data[,1],ar.data[,2:3],type="l")

write.csv(ar.data,"~/Shertzer_algae_rotifers.csv",quote=FALSE,row.names=FALSE)
