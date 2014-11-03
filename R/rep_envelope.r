# Functions to create upper and lower envelope



# L1 is a 2x2 matrix where columns are X,Y and rows are the start and end points	
cut.point.between.lines = function(L1, L2) {
	x1 = L1[1,1];x2 = L1[2,1]; y1=L1[1,2]; y2 = L1[2,2];
	w1 = L2[1,1];w2 = L2[2,1]; z1=L2[1,2]; z2 = L2[2,2];
	
	x = (((w1*x1*y2-w1*x2*y1-w1*x1*z2-w2*x1*y2+w2*x1*z1+w2*x2*y1+w1*x2*z2-w2*x2*z1))/
	    (w1*y2-w1*y1+w2*y1-w2*y2+x1*z1-x1*z2-x2*z1+x2*z2))
	    
	
	y = (((w2*y1*z1-w1*y1*z2+w1*y2*z2-w2*y2*z1+x1*y2*z1-x2*y1*z1-x1*y2*z2+x2*y1*z2))/
	    (w1*y2-w1*y1+w2*y1-w2*y2+x1*z1-x1*z2-x2*z1+x2*z2))
	c(x,y)
}


y.on.line = function(L1,x) {
	x1 = L1[1,1];x2 = L1[2,1]; y1=L1[1,2]; y2 = L1[2,2];
	return((x2*y1-x1*y2-y1*x+y2*x)/(x2-x1))
}
x.on.line = function(L1,y) {
	x1 = L1[1,1];x2 = L1[2,1]; y1=L1[1,2]; y2 = L1[2,2];
	return((x2*y1-x1*y2+x1*y-x2*y)/(y1-y2))
}

point.same.or.below = function(Xp,Yp,Xv,Yv,tol,also.below = TRUE) {
	same.below = rep(FALSE,NROW(Xp))
	if (length(Xv)==0) {
		return(same.below)
	}
	
	# If there are Xn points within tol range of Xo, we set their X value to that of Xo
	XpXv.diff = outer(Xp,Xv, function(x,y)(x-y))
	YpYv.diff = outer(Yp,Yv, function(x,y)(x-y))
	
	if (also.below ) {
		for (i in 1:NROW(XpXv.diff)) {
			same.below[i] = sum(XpXv.diff[i,] + tol > 0 & YpYv.diff[i,] - tol < 0) > 0
		}
	# Only same
	} else {
		for (i in 1:NROW(XpXv.diff)) {
			same.below[i] = sum(abs(XpXv.diff[i,]) < tol & abs(YpYv.diff[i,]) < tol) > 0
		}
	}
	same.below
}

		
# Let Y(X) be a monotone and convex function and X0<X1
# We know the slopes at X0 and X1, given by s0 and s1
# Convexity means s0>s1.
# We find the point (X,Y) where the following lines cross
# 	(X0,Y0) with slope s0
# 	(X1,Y1) with slope s1
# We also find the gap between the connection line between (X0,Y0) and (X1,Y1)
# and the calculated value of Y at position X

get.lower.cut.point = function(X0,X1,Y0,Y1,s0,s1) {
		#Om = (1/(s1-s0))(Y0-Y1-s0X0+s1X1)
    X = (1/(s1-s0))*
			  (Y0-Y1-s0*X0+s1*X1)			  
		
		# Y=(1/(s1-s0))(Y0*s1-Y1*s0+s1*s0*(W1-W0))			      
		Y =     (1/(s1-s0))*
			      (Y0*s1-Y1*s0+
			       s0*s1*(X1-X0))

    # The Y value on the direct line between (X0,Y0) and (X1,Y1)
   	Y.up = Y0+(X-X0) * ((Y1-Y0)/(X1-X0))
		res=(c(X,Y,Y.up-Y))
		names(res)=c("X","Y","gap")
		return(res)
}

get.lb.L.Y = function(i0,i1, LY.mat,tol= REP.CNT$TOL) {
	  # First check whether i0 and i1 are connected by a line (Ythin some tolerance level)

	  #B1.pred = B.seq[i0]+slope.seq[i0]*(L.seq[i1]-L.seq[i0])
	  
	  Y.pred = LY.mat[i0,"Y"]+LY.mat[i0,"slope"]*(LY.mat[i1,"L"]-LY.mat[i0,"L"])
		if (abs(Y.pred - LY.mat[i1,"Y"])<tol | abs(LY.mat[i0,"slope"]-LY.mat[i1,"slope"]) < tol) {
			return(NULL)
		}
		
		get.lower.cut.point(LY.mat[i0,"L"],LY.mat[i1,"L"],
												LY.mat[i0,"Y"],LY.mat[i1,"Y"],
												LY.mat[i0,"slope"],LY.mat[i1,"slope"])      												
}		
			


XYs.line = function(X,Y,slope,...) {
	a = Y - X*slope
	abline(a,slope,...)
}
# plot(1,1,xlim=c(0,2),ylim=c(0,2))
# XYs.line(1,1,-2)	


# Calculates the upper envelope of a piece-wise linear function (Xo,Yo) and a continous, concave
# piece-wise linear function (Xn,Yn). Returns a set of points that characterizes the new piece-wise
# linear function.
# Idea of the algorithm: 
# 1. Use approx.fun to get linear approximations of (Xo,Yo) and (Xn,Yn) that can be evaluated at every X level
# 2. Combine all points ((Xo,Xn),(Yo,Yn)) and sort increasingly in X (breaking ties with Y). 
# 3. Go through all points from left to right. 
#    3a. If actual type (o or n) changes chec
# ...
# There are two versions if slow.ind.change=TRUE we go step by step through each point
# if slow.ind.change = FALSE, we first identify the points where the index changes
# the second version is faster but not robust if points lie on top of each other
# if an error occurs, we run the slow version
LY.get.upper.envelope = function(Xo,Yo,Xn,Yn,ao = 1:length(Xo),an=0, do.lower.envelope = FALSE,
                                 collab = c("X","Y","a"),tol = REP.CNT$TOL.LY.get.upper.envelope, plot.main = "",
                                 slow.ind.change = TRUE, do.plot = REP.CNT$PLOT.DURING.SOLVE) {		
                                 
	
	#restore.point("LY.get.upper.envelope")
	
	#if (runif(1) > 0.7) {
  	#stop()
  #}
  sig = 1		
  if (do.lower.envelope) {
	  sig = -1
	  Yo = -Yo; Yn = -Yn
 	}
	names(Xo)=names(Yo)=names(Xn)=names(Yn)=NULL
 	
	# Add a last point to Xn
	if (max(Xo)>max(Xn)) {
		Xn = c(Xn,max(Xo)+1)
		Yn = c(Yn,Yn[NROW(Yn)])
	} else {
		Xo = c(Xo,max(Xn)+1)
		Yo = c(Yo,Yo[NROW(Yo)])
		# I use the overwhelming inelegant way of assigning below, because there was a very strange error when using
		# ao = c(ao,ao[NROW(ao)])
		ao = c(ao,ao[NROW(ao)])
		if (is.na(ao[NROW(ao)])) {
			warning("Envelope is.na(ao[NROW(ao)])")
			ao[NROW(ao)]=0;
		}
	}	
	
	# Add first points
	if (min(Xn)<min(Xo) & min(Yn) < min(Yo)) {
  	Xo = c(Xo[1],Xo)
  	Yo = c(min(Yn)-1,Yo)
  	ao = c(-1,ao)
	} else if (min(Xo)<min(Xn) & min(Yo) < min(Yn)) {
  	Xn = c(Xn[1],Xn)
  	Yn = c(min(Yo)-1,Yn)
	}
	
	# If there are Xn points within tol range of Xo, we set their X value to that of Xo
	XnXo.diff = outer(Xn,Xo, function(x,y)(x-y))
	reset = apply(XnXo.diff,1,function(x) {
		if (sum(x < 0 & abs(x)<tol)>0) {
	  	return(which(x < 0 & abs(x)<tol)[1])
  	} else {
	  	return(NA)
  	}
 	}) 
	Xn[!is.na(reset)] = Xo[reset[!is.na(reset)]]
	
	# Generate approximation functions
 	if (length(Xn)>=2) { 		                                                               
 		approxfun.Yn = approxfun(Xn, Yn, method="linear",yleft=-Inf, yright=max(Yn), rule = 1, ties = "ordered")
	} else {
 		approxfun.Yn = approxfun(c(Xn,Xn+1), c(Yn,Yn), method="linear",yleft=-Inf, yright=max(Yn), rule = 1, ties = "ordered")
	}
 	if (length(Xo)>=2) { 		                                                               
	 	approxfun.Yo = approxfun(Xo, Yo, method="linear",yleft=-Inf, yright=max(Yo), rule = 1, ties = "ordered")
	} else {
 		approxfun.Yo = approxfun(c(Xo,Xo+1), c(Yo,Yo), method="linear", yleft=-Inf, yright=max(Yo), rule = 1, ties = "ordered")
	}
	dist.YnYo = function(x) {approxfun.Yo(x)-approxfun.Yn(x)}
	dist.min = max(min(Xo),min(Xn))
	
	xli = range(Xo,Xn)
	yli = range(sig*Yo,sig*Yn)
	yli[2] = yli[2]+0.2*diff(yli)

	if (do.plot) {
	  plot(Xo,sig*Yo,type="l",xlim=xli,ylim=yli,main=plot.main,col="orange",ylab="",xlab="L")
	  lines(Xn,sig*Yn,col="green",lty=2)
	  points(Xo,sig*Yo,col="orange")
	  points(Xn,sig*Yn,col="green")
  }

	
	# Combine points and sort them according to X increasing and in a second dimension to Y decreasing
	TYPE.N = 2; TYPE.O = 1;
	# Merge points
	Xm = c(Xn,Xo)
	Ym = c(Yn,Yo)
	am = c(rep(an,NROW(Xn)),ao)  # Action profile associated with the point
	tm = c(rep(2,NROW(Xn)),rep(1,NROW(Xo)))  # Type of points: 1=o 2=n
	dm = rep(FALSE,NROW(Xm))  # Stores info whether point will be deleted or not
	
	# Sort all points
	ord = order(Xm,Ym)
	Xm = Xm[ord];Ym = Ym[ord];
	am = am[ord];tm = tm[ord];
	
	
	# Initialize loop
	type.opt = tm[1] 
	ind.change = which(diff(tm)!=0)+1   # Indices where type of point changes 
	approx.Y.li = list(approxfun.Yo,approxfun.Yn)
	
	is.above.other.type = rep(NA,NROW(Xm))
	rows = (tm==1)
	is.above.other.type[rows] = approx.Y.li[[2]](Xm[rows]) < Ym[rows]
	rows = (tm==2)
	is.above.other.type[rows] = approx.Y.li[[1]](Xm[rows]) < Ym[rows]

	is.below.other.type = rep(NA,NROW(Xm))
	rows = (tm==1)
	is.below.other.type[rows] = approx.Y.li[[2]](Xm[rows]) > Ym[rows]
	rows = (tm==2)
	is.below.other.type[rows] = approx.Y.li[[1]](Xm[rows]) > Ym[rows]

	rbind(Xm,Ym,tm,is.above.other.type)

	
	X.li = list(Xo,Xn)
	X.add = Y.add = a.add = numeric(0);
	
	#row = 1
    
	row = 1
	type.opt = tm[1]
	error.row = NULL
	while(row < NROW(Xm)) {
		row = row+1
		i = row		
		type.other = 3-tm[i];
	
		#if (row==10)
		#    stop()
		# The type of the actual point is different from the type of points on the upper boundary
		if (tm[i] != type.opt) {
  		
			# The old optimal type remains optimal. Just delete the actual action profile
			if (!is.above.other.type[i]) {
				dm[i]=TRUE;
				if (do.plot)
				  points(Xm[dm],sig*Ym[dm],pch="x",col="black")
				next();
			}
			
			# Check whether I am just caught on a vertical line of Xo
			# Note that we decoded such endpoints with a=NA
			if (is.na(am[i-1]) & tm[i-1] != tm[i]) {
				dm[i]=TRUE;
				if (do.plot)
				  points(Xm[dm],sig*Ym[dm],pch="x",col="black")
				next();
			}
			
			# Otherwise, we find that the optimal point type changes
			# We have to calculate a cut-point to the left to see which previous points have to be removed
			vertical.below = TRUE
			left.ind = which(tm[1:(i-1)]==tm[i])
			if (NROW(left.ind)>0) {
				left.ind = max(left.ind)
				if (Xm[left.ind] < Xm[i]) {
					vertical.below = FALSE
				}
			}
			  						
			if (!vertical.below) {
				OK = TRUE
				Xc = findzero(dist.YnYo, lower = max(Xm[left.ind],dist.min), upper = Xm[i],result.tol=tol)
				
  			# Could not find a cut point
				if (is.na(Xc)) {
  				  x = seq(max(Xm[left.ind],dist.min),Xm[i],length=1000) 
            plot(x,dist.YnYo(x),type="l")
 	          abline(v=Xm[i])
 	             				    	
  				
  				 error.row = c(error.row,row)
  				 dm[i]=TRUE;
  				 
  				 if (do.plot) {
    				 plot(Xo,sig*Yo,type="l",xlim=c(1000,20000),ylim=yli,main=plot.main,col="orange",ylab="",xlab="L")
     	       lines(Xn,sig*Yn,col="green",lty=2)
   	         abline(v=Xm[i])
   	         points(Xo,sig*Yo,col="orange")
   	         points(Xn,sig*Yn,col="green")
   	         points(Xm[ind.change],Ym[ind.change],col="red",pch=3)
           }	         
 				   error.row = c(error.row,row)
				   break();
				 } 
			} else {
				Xc = Xm[i]
			}
			Yc = approx.Y.li[[type.other]](Xc)    # Take the Y-level from the other function
		
			# We simply add the cut-point
			X.add = c(X.add,Xc)
			Y.add = c(Y.add,Yc)
			
			if (vertical.below) {
				a.add = c(a.add,NA)
			} else {
				a.add = c(a.add,am[left.ind])
			}
			if (do.plot)
			  points(Xc,sig*Yc,col="blue")
											
			type.opt = tm[i]
			
					
		}	else   { # (tm[i] == type.opt)
		
			# Check if the actual optimal type does not remain optimal.
			# We have to calculate a cut-point to the left, add the cut-point and then remove the actual point
			if (is.below.other.type[i]) {
				
				# Check whether the point ends between a vertical line
				# In this case, we keep the point
				if (Xm[i-1] == Xm[i]) next()
						
				
				left.ind = which(tm[1:(i-1)]==tm[i])
				left.ind = max(left.ind)
				
				vertical.below = TRUE
				if (Xm[left.ind] < Xm[i]) {
					vertical.below = FALSE
				}
				
				if (! vertical.below) {
					if (dm[left.ind]) {
						Xc.left = min(X.add[X.add>=Xm[left.ind]])+tol / 2
					} else {
						Xc.left = Xm[left.ind]
					}
 					Xc = findzero(dist.YnYo, lower = max(Xc.left,dist.min), upper = Xm[i],result.tol=tol*10)
  				if (is.na(Xc)) {
    				
   				  x = seq(max(Xc.left,dist.min),Xm[i],length=1000)
 
            plot(x,dist.YnYo(x),type="l")
 	          abline(v=Xm[i])
 	             				    				 				   
 				    plot(Xo,sig*Yo,type="l",xlim=c(50000,120000),ylim=yli,main=plot.main,col="orange",ylab="",xlab="L")
 	          lines(Xn,sig*Yn,col="green",lty=2)
 	          abline(v=Xm[i])
 	          points(Xo,sig*Yo,col="orange")
 	          points(Xn,sig*Yn,col="green")
	         
  				  error.row = c(error.row,row)
 				    break();
  				} 				          
				} else {
					Xc = Xm[i]
				}

				Yc = approx.Y.li[[type.opt]](Xc)    # Take the Y-level from the actual function (is above from left)
			
				# We add the cut-point
				X.add = c(X.add,Xc)
				Y.add = c(Y.add,Yc)
				if (vertical.below) {
					a.add = c(a.add,NA)
				} else {
					ind = max(which(tm != tm[i] & Xm <= Xc))
					a.add = c(a.add,am[ind])     # Take a from the other function (is above right)
				}

				if (do.plot)
  				points(Xc,sig*Yc,col="blue")
				# Delete the actual point and the points to the right
				dm[i] = TRUE
				if (do.plot)
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
				type.opt = type.other
				
			}
		}
	}    
  	
	if (!is.null(error.row)) {
 	 	warning(paste("LY.get.upper.envelope: Error in findzero. a=",m$lab.a[an],". Return original points."))
  	
  	dev.act = dev.cur()
    win.graph()
    
    plot(Xo,sig*Yo,type="l",xlim=xli,ylim=yli,main=paste("Error in findzero a=",m$lab.a[an],"  No problem if orange >= green"),col="orange",xlab="L")
	  lines(Xn,sig*Yn,col="green",lty=2)
      
    dev.set(dev.act)
    bringToTop(which = dev.act, stay = FALSE)
    mat = cbind(Xo,Yo,ao)
	  colnames(mat)=collab	
    if (do.lower.envelope) {mat[,2] = -mat[,2]}
    #stop()
    return(mat)
	}
	
	# Merge points again
	Xa = c(Xm[!dm],X.add)
	Ya = c(Ym[!dm],Y.add)
	aa = c(am[!dm],a.add)
	ord = order(Xa,Ya)
	Xa = Xa[ord];Ya = Ya[ord];
	aa = aa[ord];
	
	if (length(Xa) > 1) {
		# Check that Ya does not decrease
		if (min(diff(Ya))< (-tol)) {
			warning("LY.get.upper.envelope: Ya was decreasing by ", min(diff(Ya)), ". Point was removed.")
			store.traceback()
			stop(paste("Stopping: LY.get.upper.envelope: Ya was decreasing by ", min(diff(Ya)), ". Point was removed."))
		} 
		
	
		# Remove points where Xa and Ya are too close to each other
		rm.ind = rep(FALSE,length(Xa))
		for (i in length(Xa):2) {
			rm.ind[i] = min(abs(Xa[1:(i-1)]-Xa[i])+abs(Ya[1:(i-1)]-Ya[i])) < tol
		}
		# Remove points to the right that have same Y-level
		if (NROW(Xa)>1) {
			if (Ya[NROW(Xa)] == Ya[NROW(Xa)-1]) {
				rm.ind[NROW(Xa)] = TRUE
			}
		}
		if (sum(rm.ind)>0) {
			Xa = Xa[!rm.ind]
			Ya = Ya[!rm.ind]
			aa = aa[!rm.ind]
		}
		
		# Set a=NA for vertical lines
		diff.ind = which(abs(diff(Xa))< tol)
		if (NROW(diff.ind) >= 1) {
			Xa[diff.ind] = Xa[diff.ind+1]
			aa[diff.ind] = NA
		}		
	}	
	points(Xa,sig*Ya,col="red",lty=2)
	lines(Xa,sig*Ya,col="red")	
		
	mat = cbind(Xa,Ya,aa)
	colnames(mat)=collab	
	
  if (do.lower.envelope) {mat[,2] = -mat[,2]}
	
#   na.ind = which(is.na(mat[,3]))
#   if (sum(abs(diff(mat[na.ind-1,2]))>tol)>0) {
# 	  stop("calculate.upper.envelope: Wrong NA in action profiles. Check out")
#  	}
#   
	return(mat)
}




































LY.get.upper.envelope.without.slow.ind.change = function(Xo,Yo,Xn,Yn,ao = 1:length(Xo),an=0, do.lower.envelope = FALSE,
                                 collab = c("X","Y","a"),tol = REP.CNT$TOL.LY.get.upper.envelope, plot.main = "",
                                 slow.ind.change = TRUE) {		
                                 
	
	#restore.point("LY.get.upper.envelope")
	
	#if (runif(1) > 0.7) {
  	#stop()
  #}
  sig = 1		
  if (do.lower.envelope) {
	  sig = -1
	  Yo = -Yo; Yn = -Yn
 	}
	names(Xo)=names(Yo)=names(Xn)=names(Yn)=NULL
 	
	# Add a last point to Xn
	if (max(Xo)>max(Xn)) {
		Xn = c(Xn,max(Xo)+1)
		Yn = c(Yn,Yn[NROW(Yn)])
	} else {
		Xo = c(Xo,max(Xn)+1)
		Yo = c(Yo,Yo[NROW(Yo)])
		# I use the overwhelming inelegant way of assigning below, because there was a very strange error when using
		# ao = c(ao,ao[NROW(ao)])
		ao = c(ao,ao[NROW(ao)])
		if (is.na(ao[NROW(ao)])) {
			warning("Envelope is.na(ao[NROW(ao)])")
			ao[NROW(ao)]=0;
		}
	}	
	
	# Add first points
	if (min(Xn)<min(Xo) & min(Yn) < min(Yo)) {
  	Xo = c(Xo[1],Xo)
  	Yo = c(min(Yn)-1,Yo)
  	ao = c(-1,ao)
	} else if (min(Xo)<min(Xn) & min(Yo) < min(Yn)) {
  	Xn = c(Xn[1],Xn)
  	Yn = c(min(Yo)-1,Yn)
	}
	
	# If there are Xn points within tol range of Xo, we set their X value to that of Xo
	XnXo.diff = outer(Xn,Xo, function(x,y)(x-y))
	reset = apply(XnXo.diff,1,function(x) {
		if (sum(x < 0 & abs(x)<tol)>0) {
	  	return(which(x < 0 & abs(x)<tol)[1])
  	} else {
	  	return(NA)
  	}
 	}) 
	Xn[!is.na(reset)] = Xo[reset[!is.na(reset)]]
	
	# Generate approximation functions
 	if (length(Xn)>=2) { 		                                                               
 		approxfun.Yn = approxfun(Xn, Yn, method="linear",yleft=-Inf, yright=max(Yn), rule = 1, ties = "ordered")
	} else {
 		approxfun.Yn = approxfun(c(Xn,Xn+1), c(Yn,Yn), method="linear",yleft=-Inf, yright=max(Yn), rule = 1, ties = "ordered")
	}
 	if (length(Xo)>=2) { 		                                                               
	 	approxfun.Yo = approxfun(Xo, Yo, method="linear",yleft=-Inf, yright=max(Yo), rule = 1, ties = "ordered")
	} else {
 		approxfun.Yo = approxfun(c(Xo,Xo+1), c(Yo,Yo), method="linear", yleft=-Inf, yright=max(Yo), rule = 1, ties = "ordered")
	}
	dist.YnYo = function(x) {approxfun.Yo(x)-approxfun.Yn(x)}
	dist.min = max(min(Xo),min(Xn))
	
	xli = range(Xo,Xn)
	yli = range(sig*Yo,sig*Yn)
	yli[2] = yli[2]+0.1*diff(yli)

	plot(Xo,sig*Yo,type="l",xlim=xli,ylim=yli,main=plot.main,col="orange",xlab="L")
	lines(Xn,sig*Yn,col="green",lty=2)
	points(Xo,sig*Yo,col="orange")
	points(Xn,sig*Yn,col="green")

	
	# Combine points and sort them according to X increasing and in a second dimension to Y decreasing
	TYPE.N = 2; TYPE.O = 1;
	# Merge points
	Xm = c(Xn,Xo)
	Ym = c(Yn,Yo)
	am = c(rep(an,NROW(Xn)),ao)  # Action profile associated with the point
	tm = c(rep(2,NROW(Xn)),rep(1,NROW(Xo)))  # Type of points: 1=o 2=n
	dm = rep(FALSE,NROW(Xm))  # Stores info whether point will be deleted or not
	
	# Sort all points
	ord = order(Xm,Ym)
	Xm = Xm[ord];Ym = Ym[ord];
	am = am[ord];tm = tm[ord];
	
	rbind(Xm,Ym,tm)
	
	# Initialize loop
	type.opt = tm[1] 
	ind.change = which(diff(tm)!=0)+1   # Indices where type of point changes 
	approx.Y.li = list(approxfun.Yo,approxfun.Yn)
	X.li = list(Xo,Xn)
	X.add = Y.add = a.add = numeric(0);
	
	#row = 1
	
	if (!slow.ind.change) {
  	error.row = NULL
  	row = 0
  	while(row < NROW(ind.change)) {
  	#while(row < NROW(Xm)) {
  		row = row+1
  		i = ind.change[row]		
  		
  		if (row < NROW(ind.change)) {
  			right.i = ind.change[row+1]-1;
  		} else {
  			right.i = NROW(Xm)
  		}
  
  		type.other = tm[i-1];
  	
  		# The type of the actual point is different from the type of points on the upper boundary
  		if (tm[i] != type.opt) {
    		
  			# Check whether the optimal type changes
  		  change.type = approx.Y.li[[type.opt]](Xm[i]) < Ym[i]
  			# The old optimal type remains optimal. Just delete the actual action profile
  			if (!change.type) {
  				dm[i:right.i]=TRUE;
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
  				next();
  			}
  			
  			# Check whether I am just caught on a vertical line of Xo
  			# Note that we decoded such endpoints with a=NA
  			# Since points are ordered increasingly in Y the condition below implies that there is
  			if (is.na(am[i-1])) {
  				dm[i:right.i]=TRUE;
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
  				next();
  			}
  			
  			# Otherwise, we find that the optimal point type changes
  			# We have to calculate a cut-point to the left to see which previous points have to be removed
  			
  			vertical.below = TRUE
  			left.ind = which(tm[1:(i-1)]==tm[i])
  			if (NROW(left.ind)>0) {
  				left.ind = max(left.ind)
  				if (Xm[left.ind] < Xm[i]) {
  					vertical.below = FALSE
  				}
  			}
  			
  						
  			if (!vertical.below) {
  				OK = TRUE
  				Xc = findzero(dist.YnYo, lower = max(Xm[left.ind],dist.min), upper = Xm[i]-tol/10,result.tol=tol)
  				#abline(v=Xm[left.ind])
  				
    			# Error: could not find a curpoint
  				if (is.na(Xc)) {
    				 error.row = c(error.row,row)
 				   
#   				   plot(Xo,sig*Yo,type="l",xlim=c(1000,20000),ylim=yli,main=plot.main,col="orange",ylab="vi or Ue",xlab="L")
#   	         lines(Xn,sig*Yn,col="green",lty=2)
#   	         abline(v=Xm[i])
#   	         points(Xo,sig*Yo,col="orange")
#   	         points(Xn,sig*Yn,col="green")
#   	         points(Xm[ind.change],Ym[ind.change],col="red",pch=3)
  	         
  				   break();
  				 } 
  			} else {
  				Xc = Xm[i]
  			}
  			Yc = approx.Y.li[[type.other]](Xc)    # Take the Y-level from the other function
  		
  			# We simply add the cut-point
  			X.add = c(X.add,Xc)
  			Y.add = c(Y.add,Yc)
  			
  			if (vertical.below) {
  				a.add = c(a.add,NA)
  			} else {
  				a.add = c(a.add,am[left.ind])
  			}
  			points(Xc,sig*Yc,col="blue")
  											
  			type.opt = tm[i]
  			
  					
  		}	else   { # (tm[i] == type.opt)
  		
  			# The actual optimal type does not remain optimal.
  			# We have to calculate a cut-point to the left, add the cut-point and then remove the actual point
  			if (Ym[i]+tol < approx.Y.li[[type.other]](Xm[i])) {
  				
  				# Check whether the point ends between a vertical line
  				# In this case, we keep the point
  				if (Xm[i-1] == Xm[i]) next()
  						
  				
  				left.ind = which(tm[1:(i-1)]==tm[i])
  				left.ind = max(left.ind)
  				
  				vertical.below = TRUE
  				if (Xm[left.ind] < Xm[i]) {
  					vertical.below = FALSE
  				}
  				
  				if (! vertical.below) {
  					if (dm[left.ind]) {
  						Xc.left = min(X.add[X.add>=Xm[left.ind]])+tol / 2
  					} else {
  						Xc.left = Xm[left.ind]
  					}
  					#Xc = uniroot(dist.YnYo, interval=c(max(Xc.left,dist.min),Xm[i]-tol/10),tol=tol/10)$root
  				
  					Xc = findzero(dist.YnYo, lower = max(Xc.left,dist.min), upper = Xm[i]-tol/10,result.tol=tol)
    				if (is.na(Xc)) {
    				  error.row = c(error.row,row)
   				    next();
    				} 				          
  				} else {
  					Xc = Xm[i]
  				}
  
  				Yc = approx.Y.li[[type.opt]](Xc)    # Take the Y-level from the actual function (is above from left)
  			
  				# We add the cut-point
  				X.add = c(X.add,Xc)
  				Y.add = c(Y.add,Yc)
  				if (vertical.below) {
  					a.add = c(a.add,NA)
  				} else {
  					ind = max(which(tm != tm[i] & Xm <= Xc))
  					a.add = c(a.add,am[ind])     # Take a from the other function (is above right)
  				}
  
  				
  				points(Xc,sig*Yc,col="blue")
  				# Delete the actual point and the points to the right
  				dm[i:right.i] = TRUE
  				
  				# Undo some deletions
  				left.i = min(which(Xm > Xc & Ym > Yc))
  				if (left.i <= i-1) {
  					dm[left.i:(i-1)] = FALSE
  					points(Xm[left.i:(i-1)],sig*Ym[left.i:(i-1)],pch="x",col="yellow")
  				}
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
  				type.opt = type.other
  				
  			}
  		}
  	}
  	
  # Slow ind.change
  } else {
    
  	row = 1
  	type.opt = tm[1]
  	while(row < NROW(Xm)) {
  	#while(row < NROW(Xm)) {
  		row = row+1
  		i = row		
 		  type.other = 3-tm[i];
  	
 		  #if (row==6)
 		    #stop()
  		# The type of the actual point is different from the type of points on the upper boundary
  		if (tm[i] != type.opt) {
    		
  			# Check whether the optimal type changes
  		  change.type = approx.Y.li[[type.opt]](Xm[i]) < Ym[i]
  			# The old optimal type remains optimal. Just delete the actual action profile
  			if (change.type==FALSE) {
  				dm[i]=TRUE;
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
  				next();
  			}
  			
  			# Check whether I am just caught on a vertical line of Xo
  			# Note that we decoded such endpoints with a=NA
  			if (is.na(am[i-1]) & tm[i-1] != tm[i]) {
  				dm[i]=TRUE;
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
  				next();
  			}
  			
  			# Otherwise, we find that the optimal point type changes
  			# We have to calculate a cut-point to the left to see which previous points have to be removed
  			vertical.below = TRUE
  			left.ind = which(tm[1:(i-1)]==tm[i])
  			if (NROW(left.ind)>0) {
  				left.ind = max(left.ind)
  				if (Xm[left.ind] < Xm[i]) {
  					vertical.below = FALSE
  				}
  			}
  			  						
  			if (!vertical.below) {
  				OK = TRUE
  				Xc = findzero(dist.YnYo, lower = max(Xm[left.ind],dist.min), upper = Xm[i]-tol/10,result.tol=tol)
  				abline(v=Xc,lty=2)
  				
    			# Could not find a cut point
  				if (is.na(Xc)) {
    				 error.row = c(error.row,row)
    				 dm[i]=TRUE;
  				   break();
  				 } 
  			} else {
  				Xc = Xm[i]
  			}
  			Yc = approx.Y.li[[type.other]](Xc)    # Take the Y-level from the other function
  		
  			# We simply add the cut-point
  			X.add = c(X.add,Xc)
  			Y.add = c(Y.add,Yc)
  			
  			if (vertical.below) {
  				a.add = c(a.add,NA)
  			} else {
  				a.add = c(a.add,am[left.ind])
  			}
  			points(Xc,sig*Yc,col="blue")
  											
  			type.opt = tm[i]
  			
  					
  		}	else   { # (tm[i] == type.opt)
  		
  			# The actual optimal type does not remain optimal.
  			# We have to calculate a cut-point to the left, add the cut-point and then remove the actual point
  			if (Ym[i] < approx.Y.li[[type.other]](Xm[i])) {
  				
  				# Check whether the point ends between a vertical line
  				# In this case, we keep the point
  				if (Xm[i-1] == Xm[i]) next()
  						
  				
  				left.ind = which(tm[1:(i-1)]==tm[i])
  				left.ind = max(left.ind)
  				
  				vertical.below = TRUE
  				if (Xm[left.ind] < Xm[i]) {
  					vertical.below = FALSE
  				}
  				
  				if (! vertical.below) {
  					if (dm[left.ind]) {
  						Xc.left = min(X.add[X.add>=Xm[left.ind]])+tol / 2
  					} else {
  						Xc.left = Xm[left.ind]
  					}
   					Xc = findzero(dist.YnYo, lower = max(Xc.left,dist.min), upper = Xm[i]-tol/10,result.tol=tol)
    				if (is.na(Xc)) {
    				  error.row = c(error.row,row)
   				    break();
    				} 				          
  				} else {
  					Xc = Xm[i]
  				}
  
  				Yc = approx.Y.li[[type.opt]](Xc)    # Take the Y-level from the actual function (is above from left)
  			
  				# We add the cut-point
  				X.add = c(X.add,Xc)
  				Y.add = c(Y.add,Yc)
  				if (vertical.below) {
  					a.add = c(a.add,NA)
  				} else {
  					ind = max(which(tm != tm[i] & Xm <= Xc))
  					a.add = c(a.add,am[ind])     # Take a from the other function (is above right)
  				}
  
  				
  				points(Xc,sig*Yc,col="blue")
  				# Delete the actual point and the points to the right
  				dm[i] = TRUE
  				points(Xm[dm],sig*Ym[dm],pch="x",col="black")
  				type.opt = type.other
  				
  			}
  		}
  	}
	}
    
  	
	if (!is.null(error.row)) {
   	# Try again using the slow.ind.change algorithm
  	if (slow.ind.change == FALSE) {
    	if (do.lower.envelope) {
      	Yo = -Yo; Yn = -Yn;
      }
      dev.act = dev.cur()
      win.graph()
      
        
      

      ret = LY.get.upper.envelope(Xo,Yo,Xn,Yn,ao,an, do.lower.envelope,
                                 collab,tol, plot.main=paste("slow.ind.change a=", m$lab.a[an]) ,
                                 slow.ind.change=TRUE)     	
      dev.set(dev.act)
      bringToTop(which = dev.act, stay = FALSE)
                                 
      return(ret)
    
  	# Return original set of points
	  } else if (slow.ind.change == TRUE) {
    	warning(paste("LY.get.upper.envelope: Error in uniroot. Return original points."))
    	
    	dev.act = dev.cur()
      win.graph()
      
      plot(Xo,sig*Yo,type="l",xlim=xli,ylim=yli,
           main=paste("Error in uniroot a=",m$lab.a[an],"  No problem if orange >= green"),col="orange",xlab="L")
  	  lines(Xn,sig*Yn,col="green",lty=2)
        
      dev.set(dev.act)
      bringToTop(which = dev.act, stay = FALSE)
      stop()
      mat = cbind(Xo,Yo,ao)
  	  colnames(mat)=collab	
      if (do.lower.envelope) {mat[,2] = -mat[,2]}
      return(mat)
    }
	}
	
	# Merge points again
	Xa = c(Xm[!dm],X.add)
	Ya = c(Ym[!dm],Y.add)
	aa = c(am[!dm],a.add)
	ord = order(Xa,Ya)
	Xa = Xa[ord];Ya = Ya[ord];
	aa = aa[ord];
	
	if (length(Xa) > 1) {
		# Check that Ya does not decrease
		if (min(diff(Ya))< (-tol)) {
			warning("LY.get.upper.envelope: Ya was decreasing by ", min(diff(Ya)), ". Point was removed.")
			store.traceback()
			stop(paste("Stopping: LY.get.upper.envelope: Ya was decreasing by ", min(diff(Ya)), ". Point was removed."))
		} 
		
	
		# Remove points where Xa and Ya are too close to each other
		rm.ind = rep(FALSE,length(Xa))
		for (i in length(Xa):2) {
			rm.ind[i] = min(abs(Xa[1:(i-1)]-Xa[i])+abs(Ya[1:(i-1)]-Ya[i])) < tol
		}
		# Remove points to the right that have same Y-level
		if (NROW(Xa)>1) {
			if (Ya[NROW(Xa)] == Ya[NROW(Xa)-1]) {
				rm.ind[NROW(Xa)] = TRUE
			}
		}
		if (sum(rm.ind)>0) {
			Xa = Xa[!rm.ind]
			Ya = Ya[!rm.ind]
			aa = aa[!rm.ind]
		}
		
		# Set a=NA for vertical lines
		diff.ind = which(abs(diff(Xa))< tol)
		if (NROW(diff.ind) >= 1) {
			Xa[diff.ind] = Xa[diff.ind+1]
			aa[diff.ind] = NA
		}		
	}	
	points(Xa,sig*Ya,col="red",lty=2)
	lines(Xa,sig*Ya,col="red")	
		
	mat = cbind(Xa,Ya,aa)
	colnames(mat)=collab	
	
  if (do.lower.envelope) {mat[,2] = -mat[,2]}
	
#   na.ind = which(is.na(mat[,3]))
#   if (sum(abs(diff(mat[na.ind-1,2]))>tol)>0) {
# 	  stop("calculate.upper.envelope: Wrong NA in action profiles. Check out")
#  	}
#   
	return(mat)
}












