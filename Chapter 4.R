################################################################
## This code contains the relevant functions for my Chapter 4 ##
################################################################

#We essentially here just use the 2D code and vary input values to explore different
 #correlation lenghts (theta values here), our Sigmas, and number of training
 # and simulation points.

library(tgp)
#Function for emulation
emulator_2d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(1,1),delta=10^-5){
  nx_j <- nrow(x_j)				# number of points in x_j 
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  E_fx <- rep(beta_0,nx)			# needed for BL update
  E_D  <- rep(beta_0,nx_j)		# needed for BL update
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  # bit that makes f(x) a smooth function!
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						# BL update
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)				# BL update	
  sdD_fx <- sqrt(diag(VarD_fx))
  
  ### return the adjusted expection and just the standard deviations ###
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}

#Function for n waves of History Matching
twodimnwave <- function(n,range,nx_j,nx,z,Impcap,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea=0.05,sigma_epsilona=0,delta=10^-5,plotting="Yes"){
  # Our range is a 2x2 matrix containing two intervals, lets draw these out
  x_1min <- range[1]
  #print(x_1min)
  x_1max <- range[3]
  #print(x_1max)
  x_2min <- range[2]
  #print(x_2min)
  x_2max <- range[4]
  #print(x_2max)
  
  # par(mfrow=c(2,n)) # can't get multiple plots for some reason??
  
  # lets make our x_j and our x
  x_j <- lhs(n=nx_j,rect=range)
  #print(x_j)
  #par(mfrow=c(1,1))
  #if(plotting=="Yes"){
  #  plot(x_j,pch=16)
  #}
  x1seq <- seq(x_1min,x_1max,len=nx)
  #print(x1seq)
  x2seq <- seq(x_2min,x_2max,len=nx)
  #print(x2seq)
  x <- as.matrix(expand.grid(x1seq,x2seq))
  #print(x)
  
  # Run model
  Lat <- 300
  Inf44 <- 40
  Out44 <- 240
  fx_j <- (Lat + x_j[,1] - x_j[,2])/(Out44 - Inf44)				# more complex 2d function
  D <- fx_j	
  
  # Results wave 1
  out1 <- emulator_2d(x_j=x_j,x=x,D=D,beta_0=beta_0a,sig2=sig2a,theta=thetaa)
  ED_fx <- matrix(out1$Expect,nrow=length(x1seq))
  sdD_fx <- matrix(out1$StanDev,nrow=length(x1seq))
  
  par(mfrow=c(1,2))
  # Wave 1 expectation
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,ED_fx,color.palette=terrain.colors,xlab="Inflows 14 days",ylab="Outflows 14 days", main="Expectation - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)},lev=seq(min(ED_fx),max(ED_fx),len=20))
  }
  # W1 Implausibility
  Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilona^2 + sigma_ea^2)  )
  
  ### plot the 2d implausibility (red high: input points are bad) ###
  levs <- c(seq(0,2,0.25),max(Imp))						# levels to colour the implausibility plot
  cols <- rev(rainbow(length(levs)-1,st=0/6,en=2/6))		# colours for the implausibility plot
  
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Impausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)})
    
    par(mfrow=c(1,1))
    ### add in the non-implausible points from x in blue (the points having I(x)<2) ###
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x,pch=4,cex=0.5)})
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x[Imp<Impcap,],pch=4,cex=0.5)})
  }
  
  ###### Now we start the other n-1 waves #######
  
  #Idea: We need a 'previous wave' input and then a 'new wave' output
  # we will let lw prefix last wave items that we need for the next wave
  # and then nw prefix any new variables / outputs
  # we begin by making the 1st wave outputs our lw inputs
  # then we will make a for loop that carrys it on
  
  # Note that x always stays the same throughout waves
  lw.x_j <- x_j
  lw.fx_j <- fx_j
  lw.Imp <- Imp
  
  #begin for loop?
  k = n - 1
  for (i in c(1:k)){
    nw.x_j <- x[lw.Imp<Impcap,]
    nw.x_jwave <- rbind(lw.x_j,nw.x_j)
    nw.fx_j <- (Lat + nw.x_j[,1] - nw.x_j[,2])/(Out44 - Inf44)
    nw.D <- c(lw.fx_j,nw.fx_j)
    nw.out <- emulator_2d(x_j=nw.x_jwave,x=x,D=nw.D,beta_0=beta_0a,sig2=sig2a,theta=thetaa)
    nw.ED_fx <- matrix(nw.out$Expect,nrow=length(x1seq))                                              
    nw.sdD_fx <- matrix(nw.out$StanDev,nrow=length(x1seq))
    
    #Exp plot
    if(plotting=="Yes"){
      #filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,cex=0.5);points(x_j,pch=16,col="yellow",cex=0.5)},
      #             lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
      #filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
      #              plot.axes = {axis(1);axis(2);points(x_j,pch=16,col="yellow",cex=0.5)},
      #             lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
      par(mfrow=c(1,2))
      filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
                     plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,cex=0.5)},
                     lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
    }
    
    
    #Imp plot
    nw.Imp <- sqrt(  (nw.ED_fx - z)^2 / ( nw.sdD_fx^2 + sigma_epsilona^2 + sigma_ea^2)  )
    if(plotting=="Yes"){
      #filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,col="blue",cex=0.5)})
      #filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,col="blue",cex=0.5)})  
    }
    # then set all of these as lw before the loop resets
    lw.x_j <- nw.x_j
    lw.fx_j <- nw.fx_j
    lw.Imp <- nw.Imp
    
    # Just for my project output plot I need this next bit
    plot.x_j <- x[lw.Imp<Impcap,]
    filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
                   plot.axes = {axis(1);axis(2);points(plot.x_j,pch=4,col="blue",cex=0.5)})
    
  }
  #print ("No. of Imp points < 0.006")
  #print (length(x[nw.Imp<0.006]))
  return(list(X=x,Expect=nw.ED_fx,StanDev=nw.sdD_fx,Imp=nw.Imp))
}


xrange <- rbind(c(0,50),c(10,200)) # this is a vector

start_time <- Sys.time()
nwaveout <- twodimnwave(2,xrange,200,50,1.2,1,0,1,c(15,15),0.05,0,plotting="Yes")
########### <- function(n,range,nx_j,nx,z,Impcap,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea,sigma_epsilona,delta=10^-5){
end_time <- Sys.time()
ellapsed <- end_time - start_time
ellapsed
#head(nwaveout)

all <- cbind(nwaveout$X,as.list(nwaveout$Expect),as.list(nwaveout$Imp))
allImp <- all[all[,4]<1,]
#head(allImp)
#plot(allImp[,1],allImp[,2])
#length(allImp[,1])


#######
#######
# For our Implausibility cut off exploration, we adjust the function slightly
 # so that we can change the colouring of our implausibility plots
# We introduce "maxcol".
twodimnwave <- function(n,range,nx_j,nx,z,Impcap=1,maxcol=2,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea=0.05,sigma_epsilona=0,delta=10^-5,plotting="Yes"){
  # Our range is a 2x2 matrix containing two intervals, lets draw these out
  x_1min <- range[1]
  x_1max <- range[3]
  x_2min <- range[2]
  x_2max <- range[4]

  x_j <- lhs(n=nx_j,rect=range)
  x1seq <- seq(x_1min,x_1max,len=nx)
  x2seq <- seq(x_2min,x_2max,len=nx)
  x <- as.matrix(expand.grid(x1seq,x2seq))

  # Run model
  Lat <- 300
  Inf44 <- 40
  Out44 <- 240
  fx_j <- (Lat + x_j[,1] - x_j[,2])/(Out44 - Inf44)				# more complex 2d function
  D <- fx_j	
  
  # Results wave 1
  out1 <- emulator_2d(x_j=x_j,x=x,D=D,beta_0=beta_0a,sig2=sig2a,theta=thetaa)
  ED_fx <- matrix(out1$Expect,nrow=length(x1seq))
  sdD_fx <- matrix(out1$StanDev,nrow=length(x1seq))
  
  par(mfrow=c(1,2))
  # Wave 1 expectation
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,ED_fx,color.palette=terrain.colors,xlab="Inflows 14 days",ylab="Outflows 14 days", main="Expectation - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)},lev=seq(min(ED_fx),max(ED_fx),len=20))
  }
  # W1 Implausibility
  Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilona^2 + sigma_ea^2)  )
  
  ### plot the 2d implausibility (red high: input points are bad) ###
  levs <- c(seq(0,maxcol,0.25),max(Imp))						# levels to colour the implausibility plot
  cols <- rev(rainbow(length(levs)-1,st=0/6,en=2/6))		# colours for the implausibility plot
  
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Impausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)})
    
    par(mfrow=c(1,1))
    ### add in the non-implausible points from x in blue (the points having I(x)<2) ###
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x,pch=4,cex=0.5)})
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x[Imp<Impcap,],pch=4,cex=0.5)})
  }
  
  ###### Now we start the other n-1 waves #######
  
  #Idea: We need a 'previous wave' input and then a 'new wave' output
  # we will let lw prefix last wave items that we need for the next wave
  # and then nw prefix any new variables / outputs
  # we begin by making the 1st wave outputs our lw inputs
  # then we will make a for loop that carrys it on
  
  # Note that x always stays the same throughout waves
  lw.x_j <- x_j
  lw.fx_j <- fx_j
  lw.Imp <- Imp
  
  #begin for loop
  k = n - 1
  for (i in c(1:k)){
    nw.x_j <- x[lw.Imp<Impcap,]
    nw.x_jwave <- rbind(lw.x_j,nw.x_j)
    nw.fx_j <- (Lat + nw.x_j[,1] - nw.x_j[,2])/(Out44 - Inf44)
    nw.D <- c(lw.fx_j,nw.fx_j)
    nw.out <- emulator_2d(x_j=nw.x_jwave,x=x,D=nw.D,beta_0=beta_0a,sig2=sig2a,theta=thetaa)
    nw.ED_fx <- matrix(nw.out$Expect,nrow=length(x1seq))                                              
    nw.sdD_fx <- matrix(nw.out$StanDev,nrow=length(x1seq))
    
    #Exp plot
    if(plotting=="Yes"){
      #filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,cex=0.5);points(x_j,pch=16,col="yellow",cex=0.5)},
      #             lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
      #filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
      #              plot.axes = {axis(1);axis(2);points(x_j,pch=16,col="yellow",cex=0.5)},
      #             lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
      par(mfrow=c(1,2))
      filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
                     plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,cex=0.5)},
                     lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
    }
    
    
    #Imp plot
    nw.Imp <- sqrt(  (nw.ED_fx - z)^2 / ( nw.sdD_fx^2 + sigma_epsilona^2 + sigma_ea^2)  )
    if(plotting=="Yes"){
      #filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,col="blue",cex=0.5)})
      #filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,col="blue",cex=0.5)})  
    }
    # then set all of these as lw before the loop resets
    lw.x_j <- nw.x_j
    lw.fx_j <- nw.fx_j
    lw.Imp <- nw.Imp
    
    # Just for my project output plot I need this next bit
    plot.x_j <- x[lw.Imp<Impcap,]
    filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
                   plot.axes = {axis(1);axis(2);points(plot.x_j,pch=4,col="blue",cex=0.5)})
    
  }
  #print ("No. of Imp points < 0.006")
  #print (length(x[nw.Imp<0.006]))
  return(list(X=x,Expect=nw.ED_fx,StanDev=nw.sdD_fx,Imp=nw.Imp))
}


#Then adding in the LR values and organising our output
xrange <- rbind(c(0,50),c(10,200)) # this is a vector
start_time <- Sys.time()
nwaveout <- twodimnwave(2,xrange,200,50,1.2,0.5,1,0,1,c(15,15),0.1,0,plotting="Yes")
########### <- function(n,range,nx_j,nx,z,Impcap,maxcol,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea,sigma_epsilona,delta=10^-5){
end_time <- Sys.time()
ellapsed <- end_time - start_time
ellapsed
all <- cbind(nwaveout$X,as.list(nwaveout$Expect),as.list(nwaveout$Imp))
allImp <- all[all[,4]<0.5,]
# Adding in LR column
LRcheck2 <- function(Lat, Inf14, Out14, Inf44, Out44){
  LR <- (Lat+Inf14-Out14)/(Out44-Inf44)
  #barplot(c(Lat,Inf14,Out14,Inf44,Out44), names.arg=c('Lat', 'Inf14', 'Out14', 'Inf44', 'Out44'))
  LR
}
LRcol <- rep(0,(length(allImp[,1])))
for (i in (1:length(allImp[,1]))){
  LRcol[i] <- LRcheck2(300,as.numeric(allImp[i,1]),as.numeric(allImp[i,2]),40,240)
}
LRcol <- as.matrix(LRcol,ncol=1,colnames=c("LR"))
allImp2d <- cbind(allImp,LRcol)
colnames(allImp2d) <- c('Inf1:14','Out1:14','Exp','Imp','LR')
#range(allImp2d[,5])



###########
###########
# For our LHS work, we adjust the n-wave formula to contain different LHS methods as follows:
library(lhs)
twodimnwave <- function(n,range,nx_j,nx,z,Impcap,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea=0.05,sigma_epsilona=0,delta=10^-5,plotting="Yes"){
  # Our range is a 2x2 matrix containing two intervals, lets draw these out
  x_1min <- range[1]
  x_1max <- range[3]
  x_2min <- range[2]
  x_2max <- range[4]
  
  # lets make our x_j and our x
  #x_j <- lhs(n=nx_j,rect=range)
  #x_j <- randomLHS(n=nx_j,2) #this is because we are in 2D
  #x_j <- maximinLHS(n=nx_j,2) #this is because we are in 2D
  #x_j <- optimumLHS(n=nx_j,2) #this is because we are in 2D
  x_j <- improvedLHS(n=nx_j,2) #this is because we are in 2D
  
  # We then need to transform our LHS sample to the right range
  x_j[,1] <- x_j[,1]*(x_1max - x_1min) + x_1min
  x_j[,2] <- x_j[,2]*(x_2max - x_2min) + x_2min
  
  
  x1seq <- seq(x_1min,x_1max,len=nx)
  x2seq <- seq(x_2min,x_2max,len=nx)
  x <- as.matrix(expand.grid(x1seq,x2seq))
  
  # Run model
  Lat <- 300
  Inf44 <- 40
  Out44 <- 240
  fx_j <- (Lat + x_j[,1] - x_j[,2])/(Out44 - Inf44)				# more complex 2d function
  D <- fx_j	
  
  # Results wave 1
  out1 <- emulator_2d(x_j=x_j,x=x,D=D,beta_0=beta_0a,sig2=sig2a,theta=thetaa)
  ED_fx <- matrix(out1$Expect,nrow=length(x1seq))
  sdD_fx <- matrix(out1$StanDev,nrow=length(x1seq))
  
  par(mfrow=c(1,2))
  # Wave 1 expectation
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,ED_fx,color.palette=terrain.colors,xlab="Inflows 14 days",ylab="Outflows 14 days", main="Expectation - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)},lev=seq(min(ED_fx),max(ED_fx),len=20))
  }
  
  # W1 Implausibility
  Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilona^2 + sigma_ea^2)  )
  
  ### plot the 2d implausibility (red high: input points are bad) ###
  levs <- c(seq(0,2,0.25),max(Imp))						# levels to colour the implausibility plot
  cols <- rev(rainbow(length(levs)-1,st=0/6,en=2/6))		# colours for the implausibility plot
  
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Impausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)})
    
    par(mfrow=c(1,1))
    ### add in the non-implausible points from x in blue (the points having I(x)<2) ###
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x,pch=4,cex=0.5)})
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x[Imp<Impcap,],pch=4,cex=0.5)})
  }
  
  ###### Now we start the other n-1 waves #######
  
  #Idea: We need a 'previous wave' input and then a 'new wave' output
  # we will let lw prefix last wave items that we need for the next wave
  # and then nw prefix any new variables / outputs
  # we begin by making the 1st wave outputs our lw inputs
  # then we will make a for loop that carrys it on
  
  # Note that x always stays the same throughout waves
  lw.x_j <- x_j
  lw.fx_j <- fx_j
  lw.Imp <- Imp
  
  #begin for loop?
  k = n - 1
  for (i in c(1:k)){
    nw.x_j <- x[lw.Imp<Impcap,]
    nw.x_jwave <- rbind(lw.x_j,nw.x_j)
    nw.fx_j <- (Lat + nw.x_j[,1] - nw.x_j[,2])/(Out44 - Inf44)
    nw.D <- c(lw.fx_j,nw.fx_j)
    nw.out <- emulator_2d(x_j=nw.x_jwave,x=x,D=nw.D,beta_0=beta_0a,sig2=sig2a,theta=thetaa)
    nw.ED_fx <- matrix(nw.out$Expect,nrow=length(x1seq))                                              
    nw.sdD_fx <- matrix(nw.out$StanDev,nrow=length(x1seq))
    
    #Exp plot
    if(plotting=="Yes"){
      #filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,cex=0.5);points(x_j,pch=16,col="yellow",cex=0.5)},
      #             lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
      #filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
      #              plot.axes = {axis(1);axis(2);points(x_j,pch=16,col="yellow",cex=0.5)},
      #             lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
      par(mfrow=c(1,2))
      filled.contour(x1seq,x2seq,nw.ED_fx,color.palette=terrain.colors,xlab="Inflows 14 Days",ylab="Outflows 14 Days",main=c("Expectation - Wave",(i+1)),
                     plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,cex=0.5)},
                     lev=seq(min(nw.ED_fx),max(nw.ED_fx),len=20))
    }
    
    
    #Imp plot
    nw.Imp <- sqrt(  (nw.ED_fx - z)^2 / ( nw.sdD_fx^2 + sigma_epsilona^2 + sigma_ea^2)  )
    if(plotting=="Yes"){
      #filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,col="blue",cex=0.5)})
      #filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
      #              plot.axes = {axis(1);axis(2);points(nw.x_j,pch=4,col="blue",cex=0.5)})  
    }
    # then set all of these as lw before the loop resets
    lw.x_j <- nw.x_j
    lw.fx_j <- nw.fx_j
    lw.Imp <- nw.Imp
    
    # Just for my project output plot I need this next bit
    plot.x_j <- x[lw.Imp<Impcap,]
    filled.contour(x1seq,x2seq,nw.Imp,col=cols,lev=levs,xlab="Inflows 14 days",ylab="Outflows 14 days",main=c("Wave",(i+1), "Implausibility"),
                   plot.axes = {axis(1);axis(2);points(plot.x_j,pch=4,col="blue",cex=0.5)})
    
  }
  #print ("No. of Imp points < 0.006")
  #print (length(x[nw.Imp<0.006]))
  return(list(X=x,Expect=nw.ED_fx,StanDev=nw.sdD_fx,Imp=nw.Imp))
}
