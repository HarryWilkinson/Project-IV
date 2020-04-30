################################################################
## This code contains the relevant functions for my Chapter 2 ##
################################################################

#Function to check any LR output
LRcheck <- function(Lat, Inf14, Out14, Inf44, Out44){
  LR <- (Lat+Inf14-Out14)/(Out44-Inf44)
  barplot(c(Lat,Inf14,Out14,Inf44,Out44), names.arg=c('Lat', 'Inf14', 'Out14', 'Inf44', 'Out44'))
  c(LR, Lat+Inf14-Out14, Out44-Inf44)
}

#Function to take Normal samples of LR outputs, fixed Lat value
nrep <- function(n, mLat, mInf14, mOut14, mInf44, mOut44, sdperc){ #taking mean values, n reps, sd percentage of mean
  Inf14 = rnorm(n, mInf14, mInf14*sdperc) 
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  LR <- (mLat + Inf14 - Out14) / (Out44 - Inf44)
  LR
}

#Function to find the required Inf14 value given the other variables values
invrepcrit <- function(LR=1, mLat, mOut14, mInf44, mOut44){
  Inf14 <- LR*(mOut44-mInf44) - mLat + mOut14
  Inf14
}

#Function to give us the required Inf14 value for some Normal inputs, fixed LR and Lat
invrep <- function(n, LR, mLat, mOut14, mInf44, mOut44, sdperc){ #taking mean values, n reps, sd percentage of mean
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14 <- LR*(Out44-Inf44) - mLat + Out14
  Inf14
}

#Function to display the required Inf14 values given a Uniform LR target
lrdist <- function(n, mLR, mLat, mOut14, mInf44, mOut44, sdperc,k){ #taking mean values, n reps, sd percentage of mean
  # K is the amount the bank can change Inf14 by as a percentage of the Inf value
  LR = runif(n, mLR-0.05*mLR, mLR+0.05*mLR )
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14req <- LR*(Out44-Inf44) - mLat + Out14 # this value will get our LR exactly
  Inf14crit <- matrix(ncol=n, nrow=1) # this is the minimum required Inf14 for LR to be >= 100%
  for (i in 1:n){
    Inf14crit[,i] <- invrepcrit(1,mLat,Out14[i],Inf44[i],Out44[i])
  }
  x <- c(1:n)
  plot(x,Inf14req,type="l",col="green",xlab="Trials",ylab="Inflows in days 1 to 14") # So green line will give us our Inf to hit LR
  lines(x,Inf14crit,col="blue") # Blue will give us the value which keeps LR >= 100
  Inf14act <- mLR*(mOut44-mInf44) - mLat + mOut14 
  #this would be the usual/average required Inf14, so we will say we can move
  # slightly up or down from this, say +- 10% of this value is possible
  abline(h=Inf14act*(1-k),col="red")
  abline(h=Inf14act*(1+k),col="red")
}

#Alternative formula for the above, counts the number of breaches and can sort our graphs
lrdist2 <- function(n, mLR, mLat, mOut14, mInf44, mOut44, sdperc,k){ #taking mean values, n reps, sd percentage of mean
  # K the percentage amount the bank can change Inf14 by
  LR = runif(n, mLR-0.05*mLR, mLR+0.05*mLR )
  #LR = rep(mLR,n)
  Out14 = rnorm(n, mOut14, mOut14*sdperc)
  Inf44 = rnorm(n, mInf44, mInf44*sdperc)
  Out44 = rnorm(n, mOut44, mOut44*sdperc)
  Inf14req <- LR*(Out44-Inf44) - mLat + Out14
  Inf14crit <- matrix(ncol=n, nrow=1)
  for (i in 1:n){
    Inf14crit[,i] <- invrepcrit(1,mLat,Out14[i],Inf44[i],Out44[i])
  }
  x <- c(1:n)
  Inf14act <- mLR*(mOut44-mInf44) - mLat + mOut14 
  # Here is the unsorted version
  plot(x,Inf14req,type="l",col="green",xlab="Trials",ylab="Inflows in days 1 to 14",cex=0.5)
  lines(x,Inf14crit,col="blue",cex=0.5)
  abline(h=Inf14act*(1-k),col="red")
  abline(h=Inf14act*(1+k),col="red")
  
  # This sorts the req values and matches the crit to their corresponding req
  test = sort(Inf14req,index.return=TRUE)
  plot(x,test$x,type="l",col="green",xlab="Trials sorted by Inflows Required Values",ylab="Inflows in days 1 to 14",cex=0.5) # So green line will give us our Inf to hit LR
  lines(x,Inf14crit[test$ix],col="blue",cex=0.5) # Blue will give us the value which takes LR <100
  
  #The below is unsorted 
  #plot(x,Inf14req,type="l",col="green",xlab="Trials",ylab="Inflows in days 1 to 14") # So green line will give us our Inf to hit LR
  #lines(x,Inf14crit,col="blue") # Blue will give us the value which takes LR <100
  #   
  Inf14act <- mLR*(mOut44-mInf44) - mLat + mOut14 
  #this would be the usual/average required Inf14, so will say we can move
  # slightly up or down from this, say +- 10% of this value is possible
  abline(h=Inf14act*(1-k),col="red")
  abline(h=Inf14act*(1+k),col="red")
  
  greenv <- matrix(ncol=n,nrow=1)
  bluev <- matrix(ncol=n,nrow=1)
  for (i in 1:n){
    if (Inf14req[i] > Inf14act*(1+k)){
      greenv[,i] <- 1
    } else {
      greenv[,i] <- 0
    }
  }
  for (j in 1:n){
    if (Inf14crit[,j] > Inf14act*(1+k)){
      bluev[,j] <- 1
    } else {
      bluev[,j] <- 0
    }
  }
  totg <- sum(greenv) #this is how many times we cannot hit LR = 120%
  totb <- sum(bluev) #this is how many times we cannot hit LR = 100%
  print(totb)
  print(totg)
  return(c(totb,totg))
  #plot(Inf14req,Inf14crit)
}

#Function to draw a bivariate pair of random variables
bvn <- function(n,mu1,mu2,sd1,sd2,r){
  mu <- c(mu1,mu2)
  sigma <- matrix(c(sd1^2, sd1*sd2*r, sd1*sd2*r, sd2^2), nrow=2)
  A <- t(chol(sigma))
  Z <- matrix(rnorm(2*n),2,n) #just N(1,0)
  bvnout <- t(A %*% Z) + matrix(rep(mu,n), byrow=TRUE,ncol=2)
  bvnout
}

#Function to draw n pairs of bivariate random variables
bvnrep <- function(n, mLat, mInf14, mOut14, mInf44, mOut44, Infr, Outr, sdperc){ #taking mean values, n reps, sd percentage of mean
  #and correlations Infr, Outf 
  bvninf <- bvn(n,mInf14,mInf44,mInf14*sdperc,mInf44*sdperc,Infr) #first col is Inf14, second Inf44
  bvnout <- bvn(n,mOut14,mOut44,mOut14*sdperc,mOut44*sdperc,Outr) # similarly
  #print(bvninf)
  #print(bvnout)
  LR <- matrix(ncol=n, nrow=1)
  for (i in 1:n){
    LR[,i] <- (mLat + bvninf[i,1] - bvnout[i,1])/ (bvnout[i,2] - bvninf[i,2])
  }
  LR
}

#Function for non-linear SDs, based upon Lat value
nrep2 <- function(n, mLat, mInf14, mOut14, mInf44, mOut44, sdperc1, sdperc2){ #taking mean values, n reps, sd percentage of mean
  Lat <- runif(n,mLat-mLat*0.3,mLat+mLat*0.3)
  Inf14 <- rep(0,n)
  Out14 <- rep(0,n)
  Inf44 <- rep(0,n)
  Out44 <- rep(0,n)
  
  for (i in 1:n){
    if (Lat[i] < mLat){
      Inf14[i] = rnorm(1, mInf14, mInf14*sdperc1) 
      Out14[i] = rnorm(1, mOut14, mOut14*sdperc1)
      Inf44[i] = rnorm(1, mInf44, mInf44*sdperc1)
      Out44[i] = rnorm(1, mOut44, mOut44*sdperc1)
    }
    
    if (Lat[i] > mLat){
      Inf14[i] = rnorm(1, mInf14, mInf14*sdperc2) 
      Out14[i] = rnorm(1, mOut14, mOut14*sdperc2)
      Inf44[i] = rnorm(1, mInf44, mInf44*sdperc2)
      Out44[i] = rnorm(1, mOut44, mOut44*sdperc2)
    }
  }
  
  LR <- (Lat + Inf14 - Out14) / (Out44 - Inf44)
  out <- cbind(LR,Lat)
  out
  #return(matrix(c(LR,Lat),nrow=n,ncol=2))
}

#Function for drawing instead Inf14 from a Gamma distribution
gamrepInf14 <- function(n, mLat, mInf14a, mInf14b, mOut14, mInf44, mOut44){ #taking mean values, n reps, sd percentage of mean
  # We now keep Lat as constant, and take Alpha & Beta parameters for our input Gammas
  Inf14 = rgamma(n,shape=mInf14a,rate=mInf14b)
  LR <- (mLat + Inf14 - mOut14) / (mOut44 - mInf44)
  LR
}

#Function for drawing n LR samples where all of our variables from a Gamma, with Lat fixed
gamrep <- function(n, beta, mLat, mInf14a, mOut14a, mInf44a, mOut44a){ #taking mean values, n reps, sd percentage of mean
  # We now keep Lat as constant, and take Beta the same for all input Gammas
  # then we take the alpha values so that our mean LR is 1.2
  Inf14 = rgamma(n,shape=mInf14a,rate=beta)
  Out14 = rgamma(n,shape=mOut14a,rate=beta)
  Inf44 = rgamma(n,shape=mInf44a,rate=beta)
  Out44 = rgamma(n,shape=mOut44a,rate=beta)
  LR <- (mLat + Inf14 - Out14) / (Out44 - Inf44)
  LR
}



