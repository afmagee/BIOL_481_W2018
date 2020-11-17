library(MASS)
library(RColorBrewer)
library(car)

# cull data to include only breeds.of.interest
cull.data <- function(dat, breeds.of.interest) {
  dat <- dat[dat$breed %in% breeds.of.interest,]
  dat <- droplevels(dat)
  return(dat)
}


# density graphic for one trait
one.trait <- function(dat, breeds.of.interest, trait) {
  # recover()
  dat <- cull.data(dat, breeds.of.interest)
  
  # The column(s) in the data frame we want
  trait_col <- which(names(dat) == trait)
  
  # kernel density estimation per dog breed
  k <- vector("list",length(breeds.of.interest))
  max_d <- -10 # Find the maximum height we'll plot
  which_first <- 0 # Find which has that maximum
  max_x <- -10 # upper plot limit
  min_x <- 10 # lower plot limit
  for (i in 1:length(breeds.of.interest)) {
    k[[i]] <- density(dat[,trait_col][dat$breed %in% breeds.of.interest[i]])
    if ( max(k[[i]]$y) > max_d ) {
      max_d <- max(k[[i]]$y)
      which_first <- i
    }
    if ( max(k[[i]]$x) > max_x ) {
      max_x <- max(k[[i]]$x)
    }
    if ( min(k[[i]]$x) < min_x ) {
      min_x <- min(k[[i]]$x)
    }
  }
  
  plot_order <- c(1:length(breeds.of.interest))[-which_first]
  
  # Plot first
  plot(k[[which_first]],col=brewer.pal(8,"Dark2")[which_first],xlim=c(min_x,max_x),main="",xlab=trait)
  polygon(k[[which_first]],col=paste0(brewer.pal(8,"Pastel2")[which_first],"99"),border=NA)
  
  # Plot rest over that
  for (i in plot_order) {
    lines(k[[i]],col=brewer.pal(8,"Dark2")[i])
    polygon(k[[i]],col=paste0(brewer.pal(8,"Pastel2")[i],"99"),border=NA)
  }
  legend("topleft",legend=breeds.of.interest,fill=paste0(brewer.pal(8,"Pastel2")[1:length(breeds.of.interest)],"99"),border=NA,bty="n")
  
}


# scatterplot graphic for two traits
two.traits <- function(dat, breeds.of.interest, trait1, trait2) {
  # recover()
  
  dat <- cull.data(dat, breeds.of.interest)
  
  # The column(s) in the data frame we want
  trait_col1 <- which(names(dat) == trait1)
  trait_col2 <- which(names(dat) == trait2)
  
  breed_colors <- brewer.pal(8,"Dark2")[1:length(breeds.of.interest)]
  breed_colors <- breed_colors[as.numeric(dat$breed)]
  
  plot(x=dat[,trait_col1], y=dat[,trait_col2], col=breed_colors, xlab=trait1, ylab=trait2, main="", pch=16)

  breed_colors <- paste0(brewer.pal(8,"Dark2")[1:length(breeds.of.interest)],"50")

  for (i in 1:length(breeds.of.interest)) {
    el <- car::dataEllipse(x=dat[,trait_col1][dat$breed %in% breeds.of.interest[i]], y=dat[,trait_col2][dat$breed %in% breeds.of.interest[i]],center.pch=NA,levels=0.8,draw=FALSE)
    polygon(el,col=breed_colors[i],border=NA)
  }
  
}


# graph an LDA & PCA
LDAPCA <- function(dat, breeds.of.interest) {
  par(mfrow=c(2,1),mai=c(0.85,0.85,0.05,0.05))
  LDA <- LDAgraphic(dat, breeds.of.interest)
  PCA <- PCAgraphic(dat, breeds.of.interest)
  par(mfrow=c(1,1))
}

# compare LDA & PCA on data to LDA & PCA on newdata
LDAPCAcompare <- function(dat, newdata, breeds.of.interest) {
  par(mfrow=c(2,2),mai=c(0.85,0.85,0.05,0.05),omi=c(0,0,0.1,0))
  LDA1 <- LDAgraphic(dat, breeds.of.interest)
  mtext("original data")
  LDA2 <- LDAgraphic(newdata, breeds.of.interest)
  mtext("your altered data")
  PCA1 <- PCAgraphic(dat, breeds.of.interest)
  PCA2 <- PCAgraphic(newdata, breeds.of.interest)
  
  par(mfrow=c(1,1))
}


# plots a graphic of an LDA
LDAgraphic <- function(dat, breeds.of.interest) {
  # recover()
  dat <- cull.data(dat, breeds.of.interest)
  lda <- lda(breed ~ ., dat)                # run LDA on dat
  fulldata <- dat                           # silly, but we could add points not used in the LDA
  plda <- predict(object=lda, fulldata=dat) # get 'predicted' values for fulldata
  
  # tidy up for the graphic
  prop.lda <- lda$svd^2/sum(lda$svd^2)

  breed_colors <- brewer.pal(8,"Dark2")[1:length(breeds.of.interest)]
  breed_colors <- breed_colors[as.numeric(dat$breed)]
  
  plot(x=plda$x[,1], y=plda$x[,2], col=breed_colors, xlab=paste0("LD1 (", floor(100*prop.lda[1]+0.5), "%)"), ylab=paste0("LD2 (", floor(100*prop.lda[2]+0.5), "%)"), main="", pch=16)
  
  breed_colors <- paste0(brewer.pal(8,"Dark2")[1:length(breeds.of.interest)],"50")
  
  for (i in 1:length(breeds.of.interest)) {
    el <- car::dataEllipse(x=plda$x[,1][dat$breed %in% breeds.of.interest[i]], y=plda$x[,2][dat$breed %in% breeds.of.interest[i]],center.pch=NA,levels=0.8,draw=FALSE)
    polygon(el,col=breed_colors[i],border=NA)
  }
  
}

# return graphic of a PCA
PCAgraphic <- function(dat, breeds.of.interest) {
  # recover()
  dat <- cull.data(dat, breeds.of.interest)
  pca <- prcomp(dat[, -1], center=T, scale.=T)   # run PCA on dat
  
  # tidy up for the graphic
  prop.pca <- pca$sdev^2/sum(pca$sdev^2)

  breed_colors <- brewer.pal(8,"Dark2")[1:length(breeds.of.interest)]
  breed_colors <- breed_colors[as.numeric(dat$breed)]
  
  plot(x=pca$x[,1], y=pca$x[,2], col=breed_colors, xlab=paste0("PCA (", floor(100*prop.pca[1]+0.5), "%)"), ylab=paste0("PCA (", floor(100*prop.pca[2]+0.5), "%)"), main="", pch=16)
  
  breed_colors <- paste0(brewer.pal(8,"Dark2")[1:length(breeds.of.interest)],"50")
  
  for (i in 1:length(breeds.of.interest)) {
    el <- car::dataEllipse(x=pca$x[,1][dat$breed %in% breeds.of.interest[i]], y=pca$x[,2][dat$breed %in% breeds.of.interest[i]],center.pch=NA,levels=0.8,draw=FALSE)
    polygon(el,col=breed_colors[i],border=NA)
  }
}
