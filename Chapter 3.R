################################################################
## This code contains the relevant functions for my Chapter 3 ##
################################################################


#########################
## 1D History Matching ##
#########################

#Function for emulation - This function is based upon that in the Durham University Workshop on Bayesian Uncertainty Analysis of Complex Models
emulator_1d <- function(x_j,x,D,beta_0=0,sig2=1,theta=15){ # Where sig2 is prior variance, 
  #theta is correlation length ie how we think smooth f(x) is 
  
  ### get numbers of points from supplied x_j and x (which are now arguments of the function)
  nx_j <- length(x_j)				# number of points in x_j (7 to start with)
  nx <- length(x)					# number of points we want to evaluate the emulator for (apply the BL update for) 
  x_both <- c(x,x_j)				# join both x and x_j together (just to make covariance calculation easier)
  n <- nx + nx_j					# total number of old and new points
  
  ### Define prior expectation of f(x) and the vector D = f(x_j) (i.e. before we run the function f(x) at x_j)
  E_fx <- rep(beta_0,nx)			
  E_D <- rep(beta_0,nx_j)			
  
  ### Define prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 * exp( - as.matrix(dist(x_both))^2/theta^2 )   # this is the weird bit that makes f(x) a smooth function!
  print (Sigma)
  # So now we should have Sigma a variance matrix, which we cut into three pieces for use
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]							# see slides 15, 16 and 17 of UNICAMP_Workshop_Lect2.pdf
  
  ### The Bayes Linear Update equations:
  ### They give the 'adjusted' expectation and 'adjusted' variance of f(x) given the run data vector D = f(x_j)
  ### (Note: to multiply matrices in R use %*% and to invert a matrix use solve().)
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)							# see slides 15 UNICAMP_Workshop_Lect2.pdf
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)					# see slides 15 UNICAMP_Workshop_Lect2.pdf
  
  ### Extract the diagonal from full variance matrix and square root as we just want the individual sd's.
  sdD_fx <- sqrt(diag(VarD_fx))
  
  ### return the adjusted expection and just the standard deviations ###
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}


#Function to plot our 1D results - again in reference to the DU workshop
plot_emulator_1d <- function(ED_fx,sdD_fx,z=1.201,sigma_e=0.01,sigma_epsilon=0,yloc=1.0,plottype="emul_imp"){
  
  ### just the emulator plot ###
  if(plottype=="emul"){
    ###Do emulation plot as before ###
    plot(x,ED_fx,ty="l",col=4,ylim=c(1.0,1.3),lwd=1,xlab="Input Parameter x",ylab="Emulator for f(x)")
    lines(x,ED_fx+2*sdD_fx,col=2,lwd=1)
    lines(x,ED_fx-2*sdD_fx,col=2,lwd=1)
    points(x_j,D,pch=16,cex=1.5)				# plot the 7 actual model runs D = f(x_j)
  }
  
  ### just the implausibility plot ###
  if(plottype=="imp"){
    ### Calculate the Implausibility I(x) for all the x values from above ###
    Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
    
    ###Plot the Implausibility and include the usual cutoff of 3 sigma as a green line ###
    plot(x,Imp,ty="l",ylim=c(0,30),xlab="Input Parameter x",ylab="Implausibility",lwd=1)
    abline(h=3,col=3,lwd=1)
  }
  
  ### emulator with coloured implausibility plot ###
  if(plottype=="emul_imp"){
    ###Do emulation plot as before ###
    plot(x,ED_fx,ty="l",col=4,ylim=c(1.0,1.3),lwd=1,xlab="Inflows days 1:14",ylab="Emulator for Liquidity Ratio")
    lines(x,ED_fx+2*sdD_fx,col=2,lwd=1)
    lines(x,ED_fx-2*sdD_fx,col=2,lwd=1)
    points(x_j,D,pch=16,cex=1.5)				# plot the 6 actual model runs D = f(x_j)
    
    ### Calculate the Implausibility I(x) for all the x values from above ###
    Imp <- sqrt(  (ED_fx - z)^2 / ( sdD_fx^2 + sigma_epsilon^2 + sigma_e^2)  )
    
    ### Now add the observed data z and +- 3 sigma error bars as horizontal lines ###
    abline(h=c(z,z+3*sigma_e,z-3*sigma_e),lty=c(1,2,2),lwd=0.7)
    
    ### Now add the implausibility as coloured points (at height yloc), 
    ### with red, yellow and green representing high, borderline and low implausibility ###
    points(x[Imp>4],rep(yloc,sum(Imp>4)),pch=16,col="red",cex=1.5)				# there are more advanced ways to do this!
    points(x[(1<Imp) & (Imp<4)],rep(yloc,sum((1<Imp) & (Imp<4))),pch=16,col="yellow",cex=1.5)
    points(x[Imp<1],rep(yloc,sum(Imp<1)),pch=16,col="green",cex=1.5)
  }
  print (range(x[Imp<1]))
}

#To run we then need to:
# Define points #
#x_j <- c(5,10,15,20,25,30,35)				# points where we will run our computer model f(x_j)
#nx <- 100								# number of points we want to evaluate the emulator for (apply the BL update for) 
#x <- seq(5,35,len=nx)					# the actual new input points x where we want to evaluate the emulator
# run the model #
#fx_j <- (Lat + x_j - Out14) / (Out44 - Inf44) 
#D <- fx_j								
# run and plot the emulator #
#out <- emulator_1d(x_j=x_j,x=x,D=D,beta_0=0,sig2=1,theta=15)
#plot_emulator_1d(ED_fx=out$Expect,sdD_fx=out$StanDev,z=1.201,sigma_e=0.01,sigma_epsilon=0,yloc=1.0,plottype="emul_imp")


#########################
## 2D History Matching ##
#########################
library(tgp)

#Function for 2D emulation - based on that in the workshop
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

#Function to complete n waves of 2D HM
twodimnwave <- function(n,range,nx_j,nx,z,Impcap,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea=0.05,sigma_epsilona=0,delta=10^-5,plotting="Yes"){
  # Our range is a 2x2 matrix containing two intervals, lets draw these out
  x_1min <- range[1]
  x_1max <- range[3]
  x_2min <- range[2]
  x_2max <- range[4]

  # lets make our x_j and our x
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
  levs <- c(seq(0,2,0.25),max(Imp))						# levels to colour the implausibility plot
  cols <- rev(rainbow(length(levs)-1,st=0/6,en=2/6))		# colours for the implausibility plot
  
  if(plotting=="Yes"){
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Impausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x_j,pch=16,cex=0.5)})
    
    par(mfrow=c(1,1))
    ### add in the non-implausible points from x in blue (the points having I(x)<Impcap) ###
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x[Imp<Impcap,],pch=4,cex=0.5)})
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
  #print ("No. of Imp points < 1")
  #print (length(x[nw.Imp<1]))
  return(list(X=x,Expect=nw.ED_fx,StanDev=nw.sdD_fx,Imp=nw.Imp))
}

# To organise the output we use
nwaveout <- twodimnwave(2,xrange,200,50,1.2,1,0,1,c(15,15),0.05,0,plotting="Yes")
all <- cbind(nwaveout$X,as.list(nwaveout$Expect),as.list(nwaveout$Imp))
allImp <- all[all[,4]<1,]
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




#########################
## 5D History Matching ##
#########################
library(gtools)
library(MASS)
library(tgp)	

#Function for 5D emulation
emulator_5d <- function(x_j,x,D,beta_0=0,sig2=1,theta=c(5,5,5,5,5),delta=10^-5){
  
  ### get numbers of points from supplied x_j and x (which are now arguments of the function)
  nx_j <- nrow(x_j)				# number of points in x_j 
  #print(nx_j)
  nx <- nrow(x)					# number of points we want to evaluate the emulator for (apply the BL update for)
  

  x_both <- rbind(x,x_j)			# join both x and x_j together (just to make covariance calculation easier)
  #print(head(x_both))
  n <- nx + nx_j					# total number of old and new points
  
  ### Define prior expectation of f(x) and the vector D = f(x_j) (i.e. before we run the function f(x) at x_j)
  E_fx <- rep(beta_0,nx)			# needed for BL update
  #print(head(E_fx))
  E_D  <- rep(beta_0,nx_j)		# needed for BL update
  
  ### Define prior variances and covariances of f(x) and D = f(x_j)
  Sigma <- sig2 *( (1-delta)* exp( - as.matrix(dist( t( t(x_both)/theta) )^2 ) ) + diag(delta,n))  # bit that makes f(x) a smooth function!
  # So now we should have Sigma a variance matrix, which we cut into three pieces for use in BL update
  Var_fx 	 <- Sigma[1:nx,1:nx]
  Var_D 	 <- Sigma[(nx+1):n,(nx+1):n]
  Cov_fx_D <- Sigma[1:nx,(nx+1):n]

  ### The Bayes Linear Update equations.
  ### They give the 'adjusted' expectation and 'adjusted' variance of f(x) given the run data vector D = f(x_j)
  ### (Note: to multiply matrices in R use %*% and to invert a matrix use solve().)
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)						# BL update
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)				# BL update	
  
  ### Extract the diagonal from full variance matrix and square root as we just want the individual sd's.
  sdD_fx <- sqrt(diag(VarD_fx))
  
  ### return the adjusted expection and just the standard deviations ###
  return(list(Expect=ED_fx,StanDev=sdD_fx))
}

#I didn't get around to making a single function for the following, but to complete
 #the first wave of HM in 5D we used the following:
out <- emulator_5d(x_j=x_j,x=x,D=D,beta_0=0,sig2=1,theta=c(15,15,15,15,15))
all <- cbind(out$Expect,out$StanDev,x)
colnames(all) <- c('Exp','SD','Lat','Inf1:14','Out1:14','Inf15:44','Out15:44')
ED_fx <- all[,1]
sdD_fx <- all[,2]
# Calculate the 5D Implausibility I(x) for all the x values
Implausibility <- function(Expected, SD, z, sigma_epsilon = 0.05, sigma_e = 0){
  Imp <- rep(0,length(Expected))
  for (i in 1:(length(Expected))){
    Imp[i] <- Imp[i] + sqrt(  (Expected[i] - z)^2 / ( SD[i]^2 + sigma_epsilon^2 + sigma_e^2)  )
  } 
  Imp
}
Imp <- Implausibility(ED_fx,sdD_fx,1.2,0.05,0)
all <- cbind(Imp,all)

LRcheck2 <- function(Lat, Inf14, Out14, Inf44, Out44){
  LR <- (Lat+Inf14-Out14)/(Out44-Inf44)
  #barplot(c(Lat,Inf14,Out14,Inf44,Out44), names.arg=c('Lat', 'Inf14', 'Out14', 'Inf44', 'Out44'))
  LR
}

LRv <- rep(0,(length(all[,1])))
for (i in (1:length(all[,1]))){
  LRv[i] <- LRv[i] + LRcheck2(all[i,4],all[i,5],all[i,6],all[i,7],all[i,8])
}

LRv <- as.matrix(LRv,ncol=1,colnames=c("LR"))
all <- cbind(all,LRv)


#Our 5D plots then came from the following few plot functions
plotanx2 <- function(all,xcol,numofsig=1){
  sdI <- sd(Imp)
  meI <- mean(Imp)
  elem <- unique(all[,xcol])
  print(elem)
  len = length(elem)
  par(mfrow=c(2,2))
  lab = labels(all[2,])
  plot(all[,xcol],all[,1],ylab="Implausibility",xlab=lab[xcol])
  
  r = rep(0,len) # red points from 1d case, highly implausible
  for (i in 1:length(elem)){
    r[i] = r[i] + length(c(all[,xcol][(all[,1]>(meI+numofsig*sdI)) & (all[,xcol] == elem[i])] ))
  }
  
  a = rep(0,len) # amber points, somewhere in the middle
  for (i in 1:length(elem)){
    a[i] = a[i] + length(c(all[,xcol][((meI-numofsig*sdI)<all[,1]) & (all[,1]<(meI+numofsig*sdI)) & (all[,xcol]==elem[i])]))
  }	
  
  g = rep(0,len) # green points, lowest implausibility 
  for (i in 1:length(elem)){
    g[i] = g[i] + length(c(all[,xcol][(all[,1]<(meI-numofsig*sdI))&(all[,xcol] ==elem[i])]))
  }
  rag <- c(r,a,g)
  ag <- c(a,g)
  ragnam <- c(elem,elem,elem)
  agnam <- c(elem,elem)
  gnam <- c(elem)
  barplot(rag,xlab = "Red, Amber, Green",density=c(rep(30,15)),angle=c(rep(0,5),rep(45,5),rep(90,5)),col="blue") #names.arg = ragnam
  barplot(ag, xlab = "Amber, Green",density=c(rep(30,10)),angle=c(rep(45,5),rep(90,5)),col="blue") #,names.arg = agnam
  barplot(g, xlab = "Green", names.arg=gnam, main=lab[xcol],density=c(rep(30,5)),angle=c(rep(90,5)),col="blue")
  cutoff <- c(meI-numofsig*sdI, meI, meI+numofsig*sdI)
  cutoff
}

plotanx3 <- function(all,xcol,numofsig=1){
  sdI <- sd(Imp)
  meI <- mean(Imp)
  elem <- unique(all[,xcol])
  print(elem)
  len = length(elem)
  lab = labels(all[2,])
  #plot(all[,xcol],all[,1],ylab="Implausibility",xlab=lab[xcol])
  
  g = rep(0,len) # green points, lowest implausibility 
  for (i in 1:length(elem)){
    g[i] = g[i] + length(c(all[,xcol][(all[,1]<(meI-numofsig*sdI))&(all[,xcol] ==elem[i])]))
  }
  gnam <- c(elem)
  bar <- barplot(g, xlab = "Green", names.arg=gnam, main=lab[xcol])
  bar
}


#I believe the below code should work also for the 5D case
# 5D as one big function
fivedim <- function(xrange, nx_j, nx, fx_j, z=1.2, sigma_e=0, sigma_epsilon=0.05,numofsig=1){
  x_j <- lhs(n=nx_j,rect=xrange)
  
  x1seq <- seq(xrange[1],xrange[2],len=nx)
  x2seq <- seq(xrange[3],xrange[4],len=nx)
  x3seq <- seq(xrange[5],xrange[6],len=nx)
  x4seq <- seq(xrange[7],xrange[8],len=nx)
  x5seq <- seq(xrange[9],xrange[10],len=nx)
  x <- setNames(expand.grid(x1seq,x2seq,x3seq,x4seq,x5seq),c('x1','x2','x3','x4','x5'))
  x <- as.matrix(x)
  
  D <- fx_j	
  beta_0=0
  sig2=1
  theta=c(15,15,15,15,15)
  delta=10^-5
  
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
  out <- (list(Expect=ED_fx,StanDev=sdD_fx))
  
  all <- cbind(out$Expect,out$StanDev,x)
  colnames(all) <- c('Exp','SD','x1','x2','x3','x4','x5')
  
  ED_fx <- all[,1]
  sdD_fx <- all[,2]
  
  Imp <- rep(0,length(ED_fx))
  for (i in 1:(length(ED_fx))){
    Imp[i] <- Imp[i] + sqrt(  (ED_fx[i] - z)^2 / ( sdD_fx[i]^2 + sigma_epsilon^2 + sigma_e^2)  )
  } 
  
  all <- cbind(Imp,all)
  plot(all[,1])
  
  LRcheck2 <- function(Lat, Inf14, Out14, Inf44, Out44){
    LR <- (Lat+Inf14-Out14)/(Out44-Inf44)
    #barplot(c(Lat,Inf14,Out14,Inf44,Out44), names.arg=c('Lat', 'Inf14', 'Out14', 'Inf44', 'Out44'))
    LR
  }
  
  LRv <- rep(0,(length(all[,1])))
  for (i in (1:length(all[,1]))){
    LRv[i] <- LRv[i] + LRcheck2(all[i,4],all[i,5],all[i,6],all[i,7],all[i,8])
  }
  
  LRv <- as.matrix(LRv,ncol=1,colnames=c("LR"))
  all <- cbind(all,LRv)
  colnames(all) <- c('Imp','Exp','SD','x1','x2','x3','x4','x5','LR')
  
  plotanx2 <- function(all,xcol,numofsig=1){
    sdI <- sd(Imp)
    meI <- mean(Imp)
    elem <- unique(all[,xcol])
    len = length(elem)
    par(mfrow=c(2,2))
    plot(all[,xcol],all[,1])
    
    r = rep(0,len) # red points from 1d case, highly implausible
    for (i in 1:length(elem)){
      r[i] = r[i] + length(c(all[,xcol][(all[,1]>(meI+numofsig*sdI)) & (all[,xcol] == elem[i])] ))
    }
    
    a = rep(0,len) # amber points, somewhere in the middle
    for (i in 1:length(elem)){
      a[i] = a[i] + length(c(all[,xcol][((meI-numofsig*sdI)<all[,1]) & (all[,1]<(meI+numofsig*sdI)) & (all[,xcol]==elem[i])]))
    }	
    
    g = rep(0,len) # green points, lowest implausibility 
    for (i in 1:length(elem)){
      g[i] = g[i] + length(c(all[,xcol][(all[,1]<(meI-numofsig*sdI))&(all[,xcol] ==elem[i])]))
    }
    rag <- c(r,a,g)
    ag <- c(a,g)
    barplot(rag)
    barplot(ag)
    barplot(g)
    cutoff <- c(meI-numofsig*sdI, meI, meI+numofsig*sdI)
    cutoff
  }
  
  plotanx2(all,4,1) # Lat
  plotanx2(all,5,1) # I14
  plotanx2(all,6,1) # O14
  plotanx2(all,7,1) # I44
  plotanx2(all,8,1) # O44
  
  # we now collect the most non-implausible points somewhere?
  sdI <- sd(Imp)
  meI <- mean(Imp)
  nip <- all[all[,1]<(meI-sdI),] # non-implausible points (nip)
  return (list(nip,all))
}

#I did create a 5D second wave function - though it has not been perfected as of 
 #the time of submission.


