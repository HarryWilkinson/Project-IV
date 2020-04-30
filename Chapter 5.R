################################################################
## This code contains the relevant functions for my Chapter 5 ##
################################################################

#We bring back the 2D n-wave functions to get our non-implausible space
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


# Function for n-waves of 2D HM
twodimnwave <- function(n,range,nx_j,nx,z,Impcap,beta_0a=0,sig2a=1,thetaa=c(2.5,2.5),sigma_ea=0.05,sigma_epsilona=0,delta=10^-5,plotting="Yes"){
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
                   plot.axes = {axis(1);axis(2);points(x[Imp<Impcap,],pch=4,cex=0.5)})
    filled.contour(x1seq,x2seq,Imp,col=cols,lev=levs,xlab="Inflows 14 Days",ylab="Outflows 14 Days", main="Implausibility - Wave 1",
                   plot.axes = {axis(1);axis(2);points(x[Imp<Impcap,],pch=4,cex=0.5)})
  }
  
  ###### Now we start the other n-1 waves #######
  lw.x_j <- x_j
  lw.fx_j <- fx_j
  lw.Imp <- Imp
  
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
  return(list(X=x,Expect=nw.ED_fx,StanDev=nw.sdD_fx,Imp=nw.Imp))
}


# Function to check LR values
LRcheck2 <- function(Lat, Inf14, Out14, Inf44, Out44){
  LR <- (Lat+Inf14-Out14)/(Out44-Inf44)
  #barplot(c(Lat,Inf14,Out14,Inf44,Out44), names.arg=c('Lat', 'Inf14', 'Out14', 'Inf44', 'Out44'))
  LR
}

# Function to find the weighted distance between two points
weuc <- function(x,y,w){
  n <- length(x)
  d <- rep(0,n)
  for (i in 1:n){
    d[i] <- w[i] * (x[i] - y[i])^2
  }
  out <- sqrt(sum(d))
  out
}

# Function to find the distance between all non-imp points and our point
allweuc <- function(val,Imppts,w){ #Imppts has cols I14, O14, Exp, Imp, LR
  n <- length(Imppts[,1])
  d <- rep(0,n)
  x <- val
  y <- Imppts[,c(1,2)]
  
  for (i in 1:n){
    k <- weuc(as.vector(x),as.numeric(y[i,]),w)
    d[i] <- k
  }
  out <- cbind(d,Imppts)
  out
}


## We then have the setup
xrange <- rbind(c(0,50),c(10,200)) 
nwaveout <- twodimnwave(2,xrange,200,50,1.2,1,0,1,c(15,15),0.05,0,plotting="Yes")
all <- cbind(nwaveout$X,as.list(nwaveout$Expect),as.list(nwaveout$Imp))
allImp <- all[all[,4]<1,]

LRcol <- rep(0,(length(allImp[,1])))
for (i in (1:length(allImp[,1]))){
  LRcol[i] <- LRcheck2(300,as.numeric(allImp[i,1]),as.numeric(allImp[i,2]),40,240)
}
LRcol <- as.matrix(LRcol,ncol=1,colnames=c("LR"))
allImp2d <- cbind(allImp,LRcol)
colnames(allImp2d) <- c('Inf1:14','Out1:14','Exp','Imp','LR')

head(allImp2d) # this is our non-implausible points 


#We can then find minimum distance points (to (5,95) here, and with weights (1,1)) using
out4 <- as.data.frame(allweuc(c(5,95),allImp2d,c(1,1)))
out4o <- as.data.frame(lapply(out4, unlist))
out4o[which.min(out4o$d),] # choosing min point



#We now move on to our random walks
#Function for one step in our random walk
oneday <- function(pos,steprangeI,steprangeO){
  p <- pos
  p[1] <- p[1] + runif(1,steprangeI[1],steprangeI[2])
  p[2] <- p[2] + runif(1,steprangeO[1],steprangeO[2])
  p
}

#Basic random walk function
walkies <- function(startpos,ndays,steprangeI,steprangeO){
  s <- startpos
  intervs <- 0
  for (i in 1:ndays){
    s <- oneday(s,steprangeI,steprangeO)
    LR <- LRcheck2(300,s[1],s[2],40,240)
    if (LR > 1.25 | LR < 1.15){
      out4 <- as.data.frame(allweuc(s,allImp2d,c(1,1)))
      out4o <- as.data.frame(lapply(out4, unlist))
      a <- out4o[which.min(out4o$d),] # choosing min point
      s[1] <- as.numeric(a[2])
      s[2] <- as.numeric(a[3])
      intervs <- intervs + 1
    }
  }
  return(c(s,intervs))
}

# As above but now we can track our movements and count adjustments
walkies2 <- function(startpos,ndays,steprangeI,steprangeO){
  s <- matrix(startpos,ncol=2,nrow=1)
  intervs <- 0
  LR <- LRcheck2(300,s[1],s[2],40,240)
  LRs <- matrix(LR,nrow=1,ncol=1)
  trackpos <- matrix(startpos,ncol=2,nrow=1)
  for (i in 1:ndays){
    s <- oneday(s,steprangeI,steprangeO)
    trackpos <- rbind(trackpos,s)
    LR <- LRcheck2(300,s[1],s[2],40,240)
    LRs <- rbind(LRs,LR)
    if (LR > 1.25 | LR < 1.15 | s[1] < 0 | s[2] < 0 | s[1] > 50 | s[2] > 200){
      out4 <- as.data.frame(allweuc(s,allImp2d,c(1,1)))
      out4o <- as.data.frame(lapply(out4, unlist))
      a <- out4o[which.min(out4o$d),] # choosing min point
      s[1] <- as.numeric(a[2])
      s[2] <- as.numeric(a[3])
      intervs <- intervs + 1
      trackpos <- rbind(trackpos,s)
      LR <- LRcheck2(300,s[1],s[2],40,240)
      LRs <- rbind(LRs,LR)
    }
  }
  trackpos <- cbind(trackpos,LRs)
  return(list(intervs,trackpos))
}

#We now only allow ourselves to adjust by a certain amount each time period
walkies3 <- function(startpos,ndays,steprangeI,steprangeO,adjlength=c(1,1)){
  s <- matrix(startpos,ncol=2,nrow=1)
  intervs <- 0
  LR <- LRcheck2(300,s[1],s[2],40,240)
  LRs <- matrix(LR,nrow=1,ncol=1)
  trackpos <- matrix(startpos,ncol=2,nrow=1)
  for (i in 1:ndays){
    s <- oneday(s,steprangeI,steprangeO)
    trackpos <- rbind(trackpos,s)
    LR <- LRcheck2(300,s[1],s[2],40,240)
    LRs <- rbind(LRs,LR)
    if (LR > 1.25 | LR < 1.15){
      if(s[1] < 0 | s[2] < 0 | s[1] > 50 | s[2] > 120){
        out4 <- as.data.frame(allweuc(s,allImp2d,c(1,1)))
        out4o <- as.data.frame(lapply(out4, unlist))
        a <- out4o[which.min(out4o$d),] # choosing min point
        s[1] <- as.numeric(a[2])
        s[2] <- as.numeric(a[3])
        intervs <- intervs + 1
        trackpos <- rbind(trackpos,s)
        LR <- LRcheck2(300,s[1],s[2],40,240)
        LRs <- rbind(LRs,LR)
      }
      out4 <- as.data.frame(allweuc(s,allImp2d,c(1,1)))
      out4o <- as.data.frame(lapply(out4, unlist))
      a <- out4o[which.min(out4o$d),] # choosing min point
      if (abs(s[1] - as.numeric(a[2])) < adjlength[1]){
        s[1] <- as.numeric(a[2])
      }
      if ((as.numeric(a[2])) > 0 & abs(s[1] - as.numeric(a[2])) > adjlength[1]){
        s[1] <- s[1] + adjlength[1]
      }
      if ((as.numeric(a[2])) < 0 & abs(s[1] - as.numeric(a[2])) > adjlength[1]){
        s[1] <- s[1] - adjlength[1]
      }
      
      if (abs(s[2] - as.numeric(a[3])) < adjlength[2]){
        s[2] <- as.numeric(a[3])      
      }
      if ((as.numeric(a[2])) > 0 & abs(s[2] - as.numeric(a[2])) > adjlength[1]){
        s[2] <- s[2] + adjlength[2]
      }
      if ((as.numeric(a[2])) < 0 & abs(s[2] - as.numeric(a[2])) > adjlength[1]){
        s[2] <- s[2] - adjlength[2]
      }    
      
      intervs <- intervs + 1
      trackpos <- rbind(trackpos,s)
      LR <- LRcheck2(300,s[1],s[2],40,240)
      LRs <- rbind(LRs,LR)
    }
  }
  trackpos <- cbind(trackpos,LRs)
  return(list(intervs,trackpos))
}


##We now consider the scenario where making adjustments cost the bank
#Function to track the associated loss with great LR values
lossfn <- function(lossfncoeff,transcost,LR,d){
  y <- lossfncoeff * (LR-1.25) + transcost*d
  y
}

# We need a function to explore how this lossfn function, the LRs and the costs interact
# so that we can pick appropriate values for our random walk scenario
explorer <- function(costpm,lossfncoeff,transcost){
  x1seq <- seq(0,50,len=51)
  x2seq <- seq(50,120,len=71) #only have 18 in this range from original 50 split
  x <- as.matrix(expand.grid(x1seq,x2seq))
  LRs <- matrix(nrow=1,ncol=1)
  for (i in 1:length(x[,1])){
    LRs <- rbind(LRs,LRcheck2(300,x[i,1],x[i,2],40,240))
  }
  LRs <- LRs[-1,]
  LRs <- as.matrix(LRs,nrow=1,ncol=1)
  x <- cbind(x,LRs)
  colnames(x) <- c('Inf1:14','Out1:14','LR')
  overlist <- matrix(nrow=1,ncol=3)
  underlist <- matrix(nrow=1,ncol=3)
  for (i in 1:length(x[,1])){
    if (x[i,3] > 1.25){
      overlist <- rbind(overlist,x[i,])
    }
    if (x[i,3]<1.15){
      underlist <- rbind(underlist,x[i,])
    }
  }
  overlist <- as.matrix(overlist[-1,],ncol=3)
  underlist <- as.matrix(underlist[-1,],ncol=3)
  
  # We now have two lists which contain the points with too high
  # and too low LR values respectively
  # We now want to explore the distances between these points and
  # our non-imp points which are in allImp2d
  
  # For each point we could find the closest point and hence the 
  # shortest distance back
  shortdist <- matrix(ncol=1,nrow=1)
  for(i in 1:length(overlist[,1])){
    out4 <- as.data.frame(allweuc(c(overlist[i,1],overlist[i,2]),allImp2d,c(1,1)))
    out4o <- as.data.frame(lapply(out4, unlist))
    a <- out4o[which.min(out4o$d),]
    shortdist <- rbind(shortdist,a[1,1])
  }
  shortdist <- as.matrix(shortdist[-1,],ncol=1)
  overlist <- cbind(overlist,shortdist)
  colnames(overlist) <- c('Inf1:14','Out1:14','LR','shortd')
  
  # Now that we have all of the over points, their LR and shortest dist to
  # the closest non-imp point, we can start exploring how our lossfn and
  # our costpm values impact these and when they would be corrected
  
  # We use
  # cost <- costpm * d
  # loss <- lossfn(lossfncoeff,transcost,LR,d)
  # and then if the loss becomes > cost to adjust, then we adjust
  
  # Ideally we would identify if they change or not on the graph
  # then add the points in either green (for unadjusted) or red (we adjusted)
  #points(x,y,col="red")
  #points(x,y,col="green")
  
  red <- matrix(ncol=4,nrow=1) # Changed from 6 to 4 in last run
  green <- matrix(ncol=4,nrow=1)
  rcl <- matrix(ncol=2,nrow=1)
  gcl <- matrix(ncol=2,nrow=1)
  for (i in 1:length(overlist[,1])){  
    d <- overlist[i,4]
    d <- d[[1]]
    LR <- overlist[i,3]
    LR <- LR[[1]]
    cost <- costpm * d
    loss <- lossfn(lossfncoeff,transcost,LR,d)
    if (loss >= cost){
      red <- rbind(red,overlist[i,])
      rcl <- rbind(rcl,c(cost,loss))
    }
    if (loss < cost){
      green <- rbind(green,overlist[i,])
      gcl <- rbind(gcl,c(cost,loss))
    }
    #if (loss >= cost){
    #  print (overlist[i,])
    #  vals <- cbind(overlist[i,],cost)
    #  vals <- cbind(vals,loss)
    #  print(vals)
    #  print(red)
    #  red <- rbind(red,vals)
    #}
    #if (loss < cost){
    #  print(overlist[i,])
    #  vals <- cbind(overlist[i,],cost)
    #  vals <- cbind(vals,loss)
    #  print(vals)
    #  print(green)
    #  green <- rbind(green,vals)
    #}
  }
  red <- as.matrix(red[-1,],ncol=6)
  green <- as.matrix(green[-1,],ncol=6)
  rcl <- as.matrix(rcl[-1,],ncol=2)
  gcl <- as.matrix(gcl[-1,],ncol=2)
  #print(head(green))
  #print(head(red))
  red <- cbind(red,rcl)
  green <- cbind(green,gcl)
  colnames(red) <- c('Inf1:14','Out1:14','LR','shortd','cost','loss')
  colnames(green) <- c('Inf1:14','Out1:14','LR','shortd','cost','loss')
  
  plot(allImp2d[,1],allImp2d[,2],xlab="Inflows 14 days",ylab="Outflows 14 day",pch=4,cex=0.4)
  points(red[,1],red[,2],col="red",cex=0.5)
  points(green[,1],green[,2],col="green",cex=0.5)
  #points(red[,1],red[,2],col='red')
  #points(green[,1],green[,2],col='green')
  
  
  return(list(overlist,green,red))
}


#Once we have chosen parameter values using the explorer function, we can again simulate
 # our random walk scenario
#Function for random walks considering associated costs
walkies4 <- function(startpos,ndays,steprangeI,steprangeO,costpm,lossfncoeff,transcost){
  s <- matrix(startpos,ncol=2,nrow=1)
  underintervs <- 0
  overintervs <- 0
  LR <- LRcheck2(300,s[1],s[2],40,240)
  LRs <- matrix(LR,nrow=1,ncol=1)
  trackpos <- matrix(startpos,ncol=2,nrow=1)
  overadj <- matrix(ncol=6,nrow=1)
  for (i in 1:ndays){
    s <- oneday(s,steprangeI,steprangeO)
    trackpos <- rbind(trackpos,s)
    LR <- LRcheck2(300,s[1],s[2],40,240)
    LRs <- rbind(LRs,LR)
    out4 <- as.data.frame(allweuc(s,allImp2d,c(1,1)))
    out4o <- as.data.frame(lapply(out4, unlist))
    a <- out4o[which.min(out4o$d),] # this is the closest point
    if (LR > 1.25){
      cost <- costpm * a[1] 
      d <- a[1]
      loss <- lossfn(lossfncoeff,transcost,LR,d)
      if (cost < loss){
        #if this is the case we want to make a movement
        olds <- s
        s[1] <- as.numeric(a[2])
        s[2] <- as.numeric(a[3])
        overintervs <- overintervs + 1
        trackpos <- rbind(trackpos,s)
        LR <- LRcheck2(300,s[1],s[2],40,240)
        LRs <- rbind(LRs,LR)
        overadjv <- c(olds,LR,d,cost,loss)
        overadj <- rbind(overadj,overadjv)
      }
    }
    if (LR < 1.10 | s[1] < 0 | s[2] < 50 | s[1] > 50 | s[2] > 200){
      s[1] <- as.numeric(a[2])
      s[2] <- as.numeric(a[3])
      underintervs <- underintervs + 1
      trackpos <- rbind(trackpos,s)
      LR <- LRcheck2(300,s[1],s[2],40,240)
      LRs <- rbind(LRs,LR)
    }
  }
  trackpos <- cbind(trackpos,LRs)
  return(list(underintervs,overintervs,trackpos,overadj))
}

