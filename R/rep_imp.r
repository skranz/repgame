
################################################################################################################
# Functions to get optimal equilibria, L(a), etc.
################################################################################################################


calculate.L.for.all.a = function(
  # Calculates the liquidity requirement for all action profiles
  # Liquidity requirements will be stored in m$L
  m, ignore.a = NULL)
{
  
  
	#restore.point("solve.all.a")

	if (k!="e")
  	stop(paste("solve.all.a so far only implement for equilibrium state k=\"e\""))
  	
	if (is.null(ignore.a))
  	ignore.a = rep(FALSE,m$nA) 
  
  rows = which(!ignore.a)
  counter = 0
  for (a in rows) {
    counter = counter+1
    if (is.na(m$L[a])) {  
      glp_std_basis(m$lp);
      ret = get.L(m,a, lp=m$lp)
      if (ret$status == 0) {
        m$L[a] = ret$L
      }
      print(paste(counter),"/",NROW(rows))
      flush.console()
    }
  }
    
  m
}

#' Solves a repeated game with voluntary transfers
#' 
#' @param m model created with init.game
#' @param symmetric if TRUE only the punishment state of player 1 is solved
#' @param ignore.ae either NULL or a BOOLEAN vector of size m$nA.
#' if ignore.ae[a] == TRUE action profile a will not be considered as an optimal 
#' equilibrium state action profile
#' @param ignore.a1 Similar to ignore.ae for punishment states of player 1 (used for symmetric games)
#' @param ignore.ai A list of vectors of action profiles that shall be ignored
#' @param cnt Control parameters for perfect monitoring games initialized with pmx.init.game
#' @param reset.L If TRUE (default) deletes all information stored in m created from a
#'                from a previous call to solve.game
#' @return returns the model m augmented by the field m$opt.mat which describes optimal action plans and corresponding joint equilibrium phase and punishment payoffs for all discount factors. See the Tutorial "Interactively Solving Repeated Games" for more information
solve.game = function(m, symmetric=m$symmetric, ignore.ae = NULL, ignore.a1 = NULL,ignore.ai = NULL, 
  reset.L = TRUE, cnt=NULL,...) {

  
	#restore.point("solve.game")
	
	if (m$sol.type=="pm" | m$sol.type=="grim" | m$sol.type == "Abreu") {
  	m = pm.solve.game(m,ignore.ae = ignore.ae,
  	                 ignore.a1 = ignore.a1,ignore.ai = ignore.ai,...)
  	return(m)
	}
	if (m$sol.type == "pmx") {
  	m = pmx.solve.game(m,cnt=cnt,...)
  	return(m)
  }  	
	
	if (reset.L) {
    m$L = rep(NA,m$nA);
    m$L[m$nash] = 0
    names(m$L) = m$lab.a
  }
	
	ret = get.max.Ue.L(m,ignore.a=ignore.ae)
	m = ret$m
	Ue.opt = ret$Ue.opt
	Ue.fun = ret$opt.approx
	#restore.point("Ue.opt")
		
	vi.opt = vi.fun = list();
	if (symmetric) {i.max = 1} else {i.max= m$n}
	for (i in 1:i.max) {
  	if (!is.null(ignore.ai)) {
      ret = get.min.vi.L(m,i = i, ignore.a = ignore.ai[[i]])
    } else {
      ret = get.min.vi.L(m,i = i, ignore.a = ignore.a1)
    }
 
		#store.objects(paste("v",i",.opt",sep=""))
		m = ret$m
		vi.opt[[i]] = ret$vi.opt
		vi.fun[[i]] = ret$opt.approx
	}
	
	opt.mat = get.opt.mat(m,Ue.opt,vi.opt)
	m$opt.mat.with.cuts = NULL
	m$opt.mat = opt.mat
	m$Ue.opt = Ue.opt
	m$vi.opt = vi.opt
	m$Ue.fun = Ue.fun
	m$vi.fun = vi.fun
	if (NROW(opt.mat) > 1) {
    m$V.fun = opt.approx = approxfun(opt.mat[,"L"], y = opt.mat[,"V"], method="linear",
                      yleft=Inf, yright=min(opt.mat[,"V"]), rule = 1, f = 0, ties = "ordered")
  } else {
    m$V.fun = opt.approx = approxfun(x = c(opt.mat[,"L"],opt.mat[,"L"]+1),
                                     y = c(opt.mat[,"V"],opt.mat[,"V"]),
                                     method="linear", ties = "ordered", rule = 1, f = 0,
                                     yleft=Inf, yright=min(opt.mat[,"V"]))
  }                                     

  m$r.fun = function(L) {
		r = m$Ue.fun(L)-m$V.fun(L) / L
		r[L==0] = Inf
		r
	}
	m$delta.fun = function(L) {1/(1+m$r.fun(L))}
		
	m$fun.of.L = list(Ue = m$Ue.fun, V=m$V.fun,delta=m$delta.fun, r=m$r.fun);
	for (i in 1:m$n) {
		if (i.max > i) {
			m$fun.of.L[[paste("v",i,sep="")]] = m$vi.fun[[i]];
		} else {
			m$fun.of.L[[paste("v",i,sep="")]] = m$vi.fun[[i.max]];
		}			
	}
	                    
  m$L.min = min(m$opt.mat[,"L"])
  m$L.max = max(m$opt.mat[,"L"])
	return(m)
	# To be continued.... need to select optimal action profiles and store everything in a matrix
}
	
# Internal function called by solve.game for games with imperfect monitoring
get.opt.mat = function(m,Ue.opt,vi.opt = NULL, sym = TRUE)
{
	
	#restore.point("get.opt.mat")

	if (sym) {i.max=1} else {i.max=m$n}
	
	# Initialize vectors that merge information from optimal matrices of all states.
	n.L = NROW(Ue.opt)
	for (i in 1:i.max) n.L = n.L + NROW(vi.opt[[i]]);
	all.L = all.k = all.a = all.Y = rep(NA,n.L);
	rows = 1:NROW(Ue.opt)
	all.L[rows] = Ue.opt[,"L"]; all.a[rows] = Ue.opt[,"a"]; all.Y[rows] = Ue.opt[,"Ue"];
	all.k[rows] = 0;
	L.start = min(Ue.opt[,"L"])
	
	for (i in 1:i.max) {
		rows = max(rows)+1:NROW(vi.opt[[i]])
		all.L[rows] = vi.opt[[i]][,"L"]; all.a[rows] = vi.opt[[i]][,"a"]; all.Y[rows] = vi.opt[[i]][,"vi"];
		all.k[rows] = i; 
		L.start = max(L.start,min(vi.opt[[i]][,"L"]))
	}
	
	# Order according to all.L and all.Y
	
	ord = order(all.L,all.Y)
	all.L = all.L[ord];all.k = all.k[ord];all.a = all.a[ord];
	
	rm.rows = all.L < L.start
	all.L = all.L[!rm.rows];all.k = all.k[!rm.rows]; all.a = all.a[!rm.rows];

	# Need to remove some elements if L.start > 0
	# Then min(L.k) can be different, in which case there do not exist optimal action structures for all L
	
	collab = c("L","delta","Ue","V",paste("v",1:m$n,sep=""),"ae",paste("a",1:m$n,sep=""),
	           "k","helper","jumped","r","UV","L.next","UV.next","opt","opt.next")
	vi.start = which(collab=="v1"); vi.cols = vi.start:(vi.start+m$n-1);
	ai.start = which(collab=="a1"); ai.cols = ai.start:(ai.start+m$n-1);	           
	
	mat = matrix(NA,NROW(all.L),NROW(collab))

	colnames(mat)=collab
	mat[,"L"] = all.L
	mat[,"k"] = all.k
	
	cbind(mat,all.a)
	#############################################################
	# Add Ue to opt.mat
	#############################################################
	
	# Add regular point of Ue to mat
	rows = (!is.na(all.a)) | all.k !=0; rows.Ue = !is.na(Ue.opt[,"a"])
	ind = findInterval(mat[rows,"L"],Ue.opt[rows.Ue,"L"])
	mat[rows,c("Ue","ae")] = Ue.opt[rows.Ue,c("Ue","a"),drop=FALSE][ind,]
	
	
	# Add helper points (points below jump point) of Ue to mat
	rows = is.na(all.a) & all.k == 0; rows.Ue = is.na(Ue.opt[,"a"])
	if (sum(rows) > 0) {
		ind = findInterval(mat[rows,"L"],Ue.opt[rows.Ue,"L"])
		mat[rows,c("Ue","ae")] = Ue.opt[rows.Ue,c("Ue","a"),drop=FALSE][ind,]
		mat[rows,"helper"] = TRUE
	}
		
	#Lue.jumped = c(FALSE,is.na(Ue.opt[1:(NROW(Ue.opt)-1),"a"]))
	#mat[mat[,"k"]==0,"jumped"] = Lue.jumped[ind[mat[,"k"]==0]]
	
	# Add the helper point at the same L level from where the jump starts

	#############################################################
	# Add all vi to opt.mat
	#############################################################
	

	
	for (i in 1:i.max) {
		# Add regular point of v1 to mat
		v.opt = vi.opt[[i]]
		ai.col = ai.cols
		
		rows = (!is.na(all.a)) | all.k !=1; rows.vi = !is.na(v.opt[,"a"])
		ind = findInterval(mat[rows,"L"],v.opt[rows.vi,"L"])
		mat[rows,c(vi.cols[i],ai.cols[i])] = v.opt[rows.vi,c("vi","a"),drop=FALSE][ind,]
	
		# Add helper points (points below jump point) of v1 to mat
		rows = is.na(all.a) & all.k == 1; rows.vi = is.na(v.opt[,"a"])
		if (sum(rows) > 0) {
			ind = findInterval(mat[rows,"L"],v.opt[rows.vi,"L"])
			mat[rows,c(vi.cols[i],ai.cols[i])] = v.opt[rows.vi,c("vi","a"),drop=FALSE][ind,]
			mat[rows,"helper"] = TRUE
		}
		
		mat[is.na(mat[,"helper"]),"helper"]=FALSE
	}
	if (sym) {
		for (i in 2:m$n) {
			mat[,ai.cols[i]] = get.symmetric.a(m,mat[,"a1"],source.i=1,dest.i=i)
		}
		mat[,vi.cols[-1]] = mat[,"v1"]
	}
			
			
	# Calculate V
	mat[,"V"] = rowSums(mat[,vi.cols,drop=FALSE])
	
	#############################################################
	# Remove points with same L. Those points with the lowest 
	#############################################################
	ord = order(mat[,"L"],mat[,"Ue"],-mat[,"V"])
	mat = mat[ord,]

	rm.rows = rep(FALSE,NROW(mat))
	rows = mat[,"helper"]==0
	rm.rows[rows] = duplicated(mat[rows,"L"],fromLast=TRUE)
	rows = mat[,"helper"]==1
	rm.rows[rows] = duplicated(mat[rows,"L"],fromLast=TRUE)
	# Make sure that there do not remain duplicates at L.start
	rows = mat[,"L"]==L.start
	rm.rows[rows] = duplicated(mat[rows,"L"],fromLast=TRUE)

	mat = mat[!rm.rows,,drop=FALSE]
	mat[1,"jumped"]=FALSE
	if (NROW(mat)>1) {
		mat[2:NROW(mat),"jumped"] = mat[1:(NROW(mat)-1),"helper"]
	}
	
	################################################################################################
	# Calculate r & delta. Afterwards mark points where delta does not fall if L gets lower
	################################################################################################
	
	# Calculate the critical discount rate
	mat[,"r"] = (mat[,"Ue"]-mat[,"V"]) / mat[,"L"]
	mat[mat[,"L"]==0,"r"] = Inf
	mat[,"delta"] = 1/(1+mat[,"r"])
	
	# The lowest discount factor of all action structures with strictly higher L
	min.delta.right = rev(cummin(c(1,rev(mat[,"delta"]))))[-1]
	
	
	# Only keep those elements that have a strictly lower critical discount factor than elements with higher L
	mat[,"opt"]= mat[,"delta"]<min.delta.right
	mat[,"UV"] = mat[,"Ue"]-mat[,"V"]

	if (NROW(mat)>1) {
  	mat[,"opt.next"]= c(mat[2:NROW(mat),"opt"],FALSE)
	  mat[,"L.next"] = c(mat[2:NROW(mat),"L"],Inf)
	  mat[,"UV.next"] = mat[c(2:NROW(mat),NROW(mat)),"UV"]
  } else {
  	mat[,"opt.next"]= c(FALSE)
	  mat[,"L.next"] = c(Inf)
	  mat[,"UV.next"] = mat[NROW(mat),"UV"]
  }		


	# Create labels
	if (!is.null(m$lab.a)) {
		lab = paste("(",m$lab.a[mat[,"ae"]],")",sep="")
		for (i in 1:i.max) {
			lab = paste(lab,",(",m$lab.a[mat[,ai.cols[i]]],")",sep="")
		}
		rownames(mat) = lab
		rownames(mat)[mat[,"helper"]==1]="---"
	}
	
		
	mat
	
	return(mat)	
}

get.vi.L.a = function(m,a,i, L.br=NULL, lp=m$lp, lp.is.for.a = FALSE, tol =  REP.CNT$TOL,max.calc=100, 
                      Y.max = Inf,LY.opt = NULL, opt.approx = NULL, do.plot=FALSE, always.calculate=FALSE,
                      L.multiplier.if.no.solution = 1.00001) {
	# Calculates vi(L|a). This function is iteratively called during solve.game for 
# models with imperfect public monitoring.
# You can also call it manually to calculate punishment payoffs manually

	
	
	#restore.point("get.vi.L.a")

	# Initialize lp such that player i has no liquidity
	if (m$lp.i == i) {
		warning(paste("get.vi.L.a: All liquidity is given to the punished player i=",i,
		              ". For higher numerical stability, one should give all liquidity to a player that is not punished"))
	}
	lambda.i = m$lambda[i]
  
	# Case that player i plays a best reply in his punishment
	# Punishment payoff is then given by c[a,i]
	if (m$g[a,i] == m$c[a,i]) {
		if (is.na(m$L[a])) {
			ret = get.L(m,a,lp=lp)	
			L = ret$L
		} else {
			L = m$L[a]
		}
		Y = m$c[a,i]
		LY.mat = matrix(c(L,Y,NA,1),1,4)
		colnames(LY.mat)=c("L","Y","slope","kink")
		res = list(br.i = TRUE,LY.mat=LY.mat, LY.bound.mat = NULL, max.gap=NULL, lp=lp, L =L)
		return(res)
	}

	# We need to know the liquidity requirement of the corresponding action profile where player i plays a best reply
	if (is.null(L.br)) {
		a.br = get.best.reply.ind(m,a,i)
		L.br = m$L[a.br]
		if (is.na(L.br)) {
			#warning(paste(sys.call(-1), 
			#	   ": you should calculate the restance of the action profile where player i best replies first"))
			warning("get.vi.L.a : you should calculate the resistance of the action profile where player i best replies first")
		  ret = get.L(m,a.br,lp=lp)
			L.br = ret$L
			lp.is.for.a = FALSE
		}
	}
	
	# If the linear program is not yet initialized for a, we need to initialize it again
	if (!lp.is.for.a) {
		make.ac.for.lp(m,a,lp) 
	}
	
	L = m$L[a]
	# If we do not yet know the liquidity requirement of a, we have to calculate it	
	if (is.na(L)) {
		change.lp.get.L(lp,m,a)
		ret = solve.glpk.lp(lp, call.lab=paste("get.vi.L.a _ get.L a=",m$lab.a[a],sep=""))
		
		# There could be action profiles that can never be implemented
		if (ret$status != 0) {
			warning(paste("get.vi.L.a: Action profile ", a , " " , m$lab.a[a], " can never be implemented."))
			res = list(br.i = FALSE,LY.mat=NULL, LY.bound.mat = NULL, max.gap=NULL, lp=lp, L =L)
			return(res)
		}
		L = ret$objval
	} 
	

	if (L.br <= L & (!always.calculate)) {
		warning(paste(sys.call(-1),"L(a) > L.br. This condition should be checked before the function is called."))
		res = list(br.i = FALSE,LY.mat=NULL, LY.bound.mat = NULL, max.gap=NULL, lp=lp, L = L)
	}

	nr = m$ny+1
	LY.mat = matrix(NA,nr,4)
	colnames(LY.mat)=c("L","Y","slope","kink")
	LY.bound.mat = matrix(NA,nr,6)
	colnames(LY.bound.mat)=c("L","Y","gap","proceed","i0","i1")
	
	# Calculate Y, i.e. v.i, for L add to LY.mat
	change.lp.get.pi.i(lp,m,a,i,L=L)
	str = paste("solve.glpk.lp from get.vi.L.a _ get.pi.i a=",m$lab.a[a]," L = ",L,sep="")
	ret = solve.glpk.lp(lp, call.lab=str,warning.if.no.solution = FALSE)
	if (ret$status != 0) {
  	change.lp.get.pi.i(lp,m,a,i,L=L*L.multiplier.if.no.solution)
	  ret = solve.glpk.lp(lp, call.lab=str)
	  if (ret$status == 0) {
  	  str = paste(str,"Solution found only after L has been multiplied by ",
  	                L.multiplier.if.no.solution)
	  } else {
  	  str = paste(str,"Solution not found even after L has been multiplied by ",
  	                L.multiplier.if.no.solution, " Assumed that player i gets maximal payments.")
  	  # Signs are reversed, so the term below denotes the maximal payments player i could receive
  	  ret$objval = (1-lambda.i) * L
	  }
	  print(str)
  	warning(str)

 	}
	
	
	Y = m$g[a,i] + ret$objval+L*lambda.i  # Note that the sign of the payments is reversed
	#Y = m$g[a,i] + ret$objval  # Note that the sign of the payments is reversed

	slope = lp.get.col.shadow.price(lp,1) + lambda.i
	LY.mat[1,] = 	c(L,Y,slope,TRUE)

	# Calculate Y, i.e. w.i, for L.br add to LY.mat
	# This an upper bound on the w.i that we want to try to achieve
	change.lp.get.pi.i(lp,m,a,i,L=L.br)
	ret = solve.glpk.lp(lp,call.lab=paste("get.vi.L.a _ get.pi.i(L.br) a=",m$lab.a[a]," L.br = ", L.br,sep=""))


	Y = m$g[a,i]+ret$objval  +L.br*lambda.i
	slope = lp.get.col.shadow.price(lp,1)  + lambda.i

	LY.mat[2,] = 	c(L.br,Y,slope,TRUE)

	# Calculate point for LY.bound.mat
	i0=1;i1=2
	point = get.lb.L.Y(i0,i1, LY.mat,tol)	
	stop.here = TRUE	
	if (!is.null(point)) {
		if (point["X"] > LY.mat[1,1] & point["X"] < LY.mat[2,1] &
		    point["Y"] < LY.mat[1,2] & point["Y"] > LY.mat[2,2]) {
			stop.here = FALSE
			row.low.LY = 1
			LY.bound.mat[row.low.LY,]=c(point,TRUE,i0,i1)
			if (LY.mat[i1,"Y"]> Y.max) {
				LY.bound.mat[row.low.LY,"proceed"]=FALSE
			}		
		}
	}
	if (stop.here) {
		res = list(LY.mat=LY.mat, LY.bound.mat = NULL, max.gap=0, lp=lp, L = L)
		res$LY.mat = res$LY.mat[!is.na(res$LY.mat[,1]),]
		return(res)
	}
	
	

	res = proceed.LY.mat(m=m,a=a,lp=lp,LY.mat=LY.mat,LY.bound.mat=LY.bound.mat,state=i,max.calc = max.calc,
	                     row.LY = 2, reduce = TRUE, Y.max = Y.max, LY.opt = LY.opt, opt.approx = opt.approx)
	
	if (do.plot) {                     
	  plot.LY.mat(res$LY.mat,res$LY.bound.mat,main="Y(L)",ylab="Y")
  }
	res = list(LY.mat=res$LY.mat, LY.bound.mat = res$LY.bound.mat, max.gap=res$max.gap,
	           lp=lp,br.i = FALSE, L = L)
	return(res)
}


# Calculates Ue(L|a). This function is iteratively called during solve.game for 
# models with imperfect public monitoring.
# You can also call it manually to calculate punishment payoffs manually		
get.Ue.L.a = function(m,a, lp=m$lp, tol =  REP.CNT$TOL,max.calc=1000, do.plot=FALSE,Y.min=-Inf, 
                      sym.y=NULL,sym.a=NULL,LY.opt = NULL, opt.approx = NULL) {

	
	#restore.point("get.Ue.L.a")
  	
			
	ret = get.L.and.B(m=m,a=a,lp=lp)
	L = ret$L		
	lp = ret$lp
	
	# For the warning messages
	warn.lab.a = m$lab.a[a]
	G.a = m$G[a]

	maybe.better.than.opt = TRUE
	if (!is.finite(ret$L)) {
		warning(paste("get.Ue.L.a: Action profile ", warn.lab.a, " cannot be implemented"))
		return(list(LY.mat=NULL, LY.bound.mat = NULL, max.gap=NULL, lp=lp, maybe.better.than.opt=FALSE, L=ret$L))
	}
	
	if (ret$L == ret$L.max) {
		LY.mat = matrix(c(ret$L,G.a-ret$B.max,NA,1),1,4)
		colnames(LY.mat)=c("L","Y","slope","kink")
		if (!is.null(opt.approx)) {
			maybe.better.than.opt = (opt.approx(LY.mat[1,"L"]) < LY.mat[1,"Y"])
		}
		res = list(LY.mat=LY.mat, LY.bound.mat = NULL, max.gap=0, lp=lp, 
		           maybe.better.than.opt = maybe.better.than.opt, L=L)
		
		return(res)
	}
	
	nr = m$ny+1
	LY.mat = matrix(NA,nr,4)
	colnames(LY.mat)=c("L","Y","slope","kink")
	LY.bound.mat = matrix(NA,nr,6)
	colnames(LY.bound.mat)=c("L","Y","gap","proceed","i0","i1")
	
	LY.mat[1,] = 	c(ret$L,G.a-ret$B.max,-ret$B.max.slope,TRUE)
	LY.mat[2,] = 	c(ret$L.max,G.a-ret$B.min,-ret$B.min.slope,TRUE)

	
	i0=1;i1=2;
	point = get.lb.L.Y(i0,i1, LY.mat,tol)	
	stop.here = TRUE	
	if (!is.null(point)) {
		if (point["X"] > LY.mat[1,1] & point["X"] < LY.mat[2,1] &
		    point["Y"] > LY.mat[1,2] & point["Y"] < LY.mat[2,2]) {
			stop.here = FALSE
			row.low.LY = 1
			LY.bound.mat[row.low.LY,]=c(point,TRUE,i0,i1)
			if (LY.mat[i1,"Y"] < Y.min) {
				LY.bound.mat[row.low.LY,"proceed"]=FALSE
			}		
		}
	}
	if (stop.here) {
		res = list(LY.mat=LY.mat, LY.bound.mat = NULL, max.gap=0, lp=lp, L=L)
		res$LY.mat = res$LY.mat[!is.na(res$LY.mat[,1]),]
		#plot.LY.mat(res$LY.mat,res$LY.bound.mat,main="U(L)",ylab="B")
		return(res)
	}


	res = proceed.LY.mat(m=m,a=a,lp=lp,LY.mat=LY.mat,LY.bound.mat=LY.bound.mat,
	                     state="e",max.calc = max.calc, row.LY = 2, reduce = TRUE, Y.min = Y.min,
	                     LY.opt = LY.opt, opt.approx = opt.approx)
	
	
	res$L = L
	if (is.null(res$LY.bound.mat)) {
		res$LY.mat = res$LY.mat[res$LY.mat[,"kink"]==1,]
	}
	if (do.plot) {
	  plot.LY.mat(res$LY.mat,res$LY.bound.mat,main=paste("Ue a =", m$lab.a[a]),ylab="Ue",xlab="L")
  }
	return(res)
}		



# Calculates the lower envelope vi(L). This function is called during solve.game for 
# models with imperfect public monitoring.
get.min.vi.L = function(m,i = 1,lp=m$lp, ignore.a, do.plot = REP.CNT$PLOT.DURING.SOLVE) {
	
  # restore.point("get.min.vi.L")

  # Make local copies
	c = m$c; G = m$G;	n = m$n; ny = m$ny;

		
	# Can never ignore Nash equilibria of the stage game
	if (is.null(ignore.a))
	  ignore.a = rep(FALSE,m$nA)
  ignore.a[m$nash] = FALSE

	
	# Initialize lp such that player i has no liquidity
	if (m$lp.i == i) {
			#remove.lp(lp)
			if (m$lp.i == m$n) {
				new.lp.i = 1
			} else {
				new.lp.i = m$n
			}
			ret = make.default.lp.for.m(m,i=new.lp.i)
			m$lp = ret$lp
			m$lp.info = ret$lp.info
			m$lp.i = ret$lp.info$i
			m$lambda = ret$lp.info$lambda
			lp = m$lp
	}

	# Check out the best Nash equilibrium
	if (NROW(m$nash)>0) {
		nash.ci = min(m$c[m$nash,i])
		}	else {
		nash.ci = Inf
  }	                               	
  
  # Specify order in which the action profiles shall be checked
  A.ord = order(!ignore.a,m$c[,i]<nash.ci,m$c[,i]==m$g[,i],-(m$C-m$G),-m$c[,i],decreasing=TRUE)
  cbind(A.ord,m$c[A.ord,i]==m$g[A.ord,i],m$c[A.ord,i],(m$C[A.ord]-m$G[A.ord]))

  max.a.check = sum(m$c[,i]<nash.ci & !ignore.a)
  max.a.check
	
	a.ind = 1
	a = A.ord[a.ind]
	
	# Find first profile a that can be implemented for some available liquidity L
	while(a.ind <= length(A.ord)) {
		ret = get.vi.L.a(m=m,a=a,i=i,lp=lp)
		if (is.null(ret$LY.mat)) {
			m$L[a] = Inf
			a.ind = a.ind +1
			a = A.ord[a.ind]
			next()
		} else {
			break()
		}
	}
			
	lp = ret$lp
		
	if (is.na(m$L[a])) {
		m$L[a] = ret$LY.mat[1,"L"]
	}	
	
	collab = c("L","vi","a")
	vi.opt = cbind(ret$LY.mat[,1],ret$LY.mat[,2],a)
	colnames(vi.opt) = collab

	Y.max = Inf
	# Add the best Nash equilibrium to vi.opt
	if (NROW(m$nash)>0 & !(m$G[vi.opt[1,"a"]]==m$C[vi.opt[1,"a"]])) {
		a = which(m$c[,i] == nash.ci & get.is.nash(m))[1]
		LY.mat = cbind(0,nash.ci,a)    # Note that for Nash equilibria indeed vi[a] = c[a,i] = g[a,i]
		colnames(LY.mat) = collab
		vi.opt = LY.get.upper.envelope(Xo=vi.opt[,1],Yo=vi.opt[,2],ao=vi.opt[,3],
		                               Xn = LY.mat[,1], Yn = LY.mat[,2], an = a, do.lower.envelope = TRUE,
		                               collab=collab)
		if(sum(!is.finite(vi.opt[,1]))>0) {
			warning("get.min.vi.L: NA's in vi.opt (adding Nash) This should not be. Debug!")
			stop()
		}
	                               		                                 
   	Y.max   = nash.ci
	}		                               
	
	a.ind.start = a.ind


	if (sum(!is.na(vi.opt[,1]))>1) {
		opt.approx = approxfun(vi.opt[,1], y = vi.opt[,2], method="linear",
                    yleft=Inf, yright=min(vi.opt[,2]), rule = 1, f = 0, ties = "ordered")
	} else {
		opt.approx = approxfun(c(vi.opt[1,1],vi.opt[1,1]+1), y = c(vi.opt[1,2],vi.opt[1,2]), method="linear",
                    yleft=Inf, yright=min(vi.opt[,2]), rule = 1, f = 0, ties = "ordered")
	}
	
	plot.main = paste("v",i, "(#A ", max.a.check, " / " ,m$nA, ")", sep="")
	if (do.plot) {
	  #plot.main = paste("v",i," #A=",m$nA,sep="")
	  plot(vi.opt[,1],vi.opt[,2],type="b",col="black",main=plot.main)
  }
	num.L.calc = 2

	L.br = NULL
	while (a.ind < m$nA) {
		
		a.ind = a.ind +1
		if (a.ind>max.a.check) {break();}
		a = A.ord[a.ind]

		# Can we rule out the profile without calculting L(a)?
		# We use the fact that m$C[a]-m$G[a] is a lower bound on the liquidity requirement
		
		
		#restore.point("loc")
		if (do.plot)
		  try(points(m$C[a]-m$G[a],m$c[a,i], col="yellow"))
		if (opt.approx(m$C[a]-m$G[a]) <= m$c[a,i])	next();
		#restore.point("loc")
		# Draw a point of the upper approximation
		if (do.plot)
  		try(legend("topright",legend=paste("a:", a.ind,"L:",num.L.calc),bg="lightgrey"))
		# Calculate a tighter upper bound without solving for L
		if (m$g[a,i] != m$c[a,i]) {
		  a.br = get.best.reply.ind(m,a,i)
			y = which.max(m$phi.mat[,a] / m$phi.mat[,a.br])[1]
			p.lb = (m$c[a,i]-m$g[a,i]) / m$phi.mat[y,a]
			vi.lb = m$c[a,i] + m$phi.mat[y,a.br]*p.lb
			try(points(m$C[a]-m$G[a],vi.lb, col="orange"))
			if (opt.approx(m$C[a]-m$G[a]) <= vi.lb)	next();
			# Not yet go to next
		} else {
			vi.lb = m$c[a,i]
		}
	

		
		made.lp = FALSE
		if (is.na(m$L[a])) {
			ret = get.L(m,a,lp=lp)
			lp = ret$lp
			m$L[a] = ret$L
			L = ret$L
			made.lp = TRUE
			num.L.calc = num.L.calc +1
		}
		
		# Has the action profile an infinite liquidity requriement?
		if (is.infinite(m$L[a])) next();

		#try(points(m$L[a],vi.lb, col="red"))

		# Maybe the profile is worse than the best that can already be achieved
		#if (opt.approx(m$L[a]) < m$c[a,i])	next();
		if (opt.approx(m$L[a]) <= vi.lb)	next();

		# If a.i is not a stage game best reply, we receive the liquidity requirement of the best-reply				
		if (m$g[a,i] != m$c[a,i]) {
		  a.br = get.best.reply.ind(m,a,i)
		  
		  # This code is inefficient. We should try to keep the old basis from get.L(m,a,lp=lp)
		  if (is.na(m$L[a.br])) {
				ret = get.L(m,a.br,lp=lp)
				m$L[a.br] = ret$L
				made.lp = FALSE  # This means get.vi has to reinitialize the problem for profile a
			}
			L.br = m$L[a.br]
			# Is L[a] >= L[a.br]? We then can neglect the profile a
			if ((m$g[a,i] != m$c[a,i]) & (m$L[a] >= L.br)) next();
		}		
		
  					
		# None of the ex-ante tricks to discard a worked, so we start calculating vi(L|a)
		# The function get.vi calculates vi(L|a) only in those parts where it can improve on the
		# existing envelope LY.opt
		
		ret = get.vi.L.a(m=m,a=a,i=i,L.br = L.br,lp=lp,lp.is.for.a = made.lp,
		                 Y.max = Y.max, LY.opt = vi.opt, opt.approx = opt.approx)
		
		if (sum(diff(ret$LY.mat[,2])> REP.CNT$TOL)>0) {
			warning("Error: vi increases with L Check out")
			restore.point("error.get.min.vi")
			#restore.point("error.get.min.vi")
			
			#copy.local.objects.to.global()
			ret = get.vi.L.a(m=m,a=a,i=i,L.br = L.br,lp=lp,lp.is.for.a = made.lp,
		                 Y.max = Y.max, LY.opt = vi.opt, opt.approx = opt.approx)
			if (sum(diff(ret$LY.mat[,2])>0)>0) {
				stop("Called again after stored global objects and still increase")
			} else {
					stop("Called again and no increase")
			}
		}
		                 		

    if (is.null(ret$LY.mat)) {
			warning(paste("get.min.vi.L: No L found for action profile ", a, " = ", m$lab.a[a]))
			print(paste("get.min.vi.L: No L found for action profile ", a, " = ", m$lab.a[a]))
			m$L[a] = Inf
			next()
		}
    
    if (!is.null(m$L[a])) {
			m$L[a] = ret$LY.mat[1,"L"]
		}
				
		vi.opt = LY.get.upper.envelope(Xo=vi.opt[,1],Yo=vi.opt[,2],ao=vi.opt[,3],
		                               Xn = ret$LY.mat[,1], Yn = ret$LY.mat[,2], an = a,
		                               do.lower.envelope = TRUE, collab=collab, plot.main=plot.main)                          
		if(sum(!is.finite(vi.opt[,1]))>0) {
			warning("get.min.vi.L: NA's in vi.opt This should not be. Debug!")
			stop()
		}
		
		if (sum(!is.na(vi.opt[,1]))>1) {
			opt.approx = approxfun(vi.opt[,1], y = vi.opt[,2], method="linear",
	                    yleft=Inf, yright=min(vi.opt[,2]), rule = 1, f = 0, ties = "ordered")
		} else {
			opt.approx = approxfun(c(vi.opt[1,1],vi.opt[1,1]+1), y = c(vi.opt[1,2],vi.opt[1,2]), method="linear",
	                    yleft=Inf, yright=min(vi.opt[,2]), rule = 1, f = 0, ties = "ordered")
		}
    
    #plot.Ue.opt(vi.opt)
 		#a.ind = a.ind +1    	                               
	}
	

	if (do.plot) {
	  plot.Ue.opt(vi.opt, main=plot.main)
  }
   
	rownames(vi.opt) = rownames(m$g)[vi.opt[,"a"]]
	return(ret=list(m=m, opt.approx = opt.approx, vi.opt = vi.opt, i = i))
}

# Function that creates the upper envelope Ue(L).
# This function is called during solve.game
get.max.Ue.L = function(m,lp=m$lp, sym.a.first = m$symmetric, ignore.a = NULL, do.plot = REP.CNT$PLOT.DURING.SOLVE) {
	
	#restore.point("get.max.Ue.L")
	
	
  # Make local copies
	c = m$c; G = m$G;	n = m$n; ny = m$ny;
	
	
	# Can never ignore Nash equilibria of the stage game
	if (is.null(ignore.a))
	  ignore.a = rep(FALSE,m$nA)
	  
  ignore.a[m$nash] = FALSE
	
	# Check out the best Nash equilibrium to Ue.opt
	if (NROW(m$nash)>0) {
		nash.G = max(m$G[m$nash])
	}	else {
		nash.G = -Inf
  }	                               	
	max.a.check = sum(m$G>nash.G & !ignore.a)
  max.a.check
  if (sym.a.first) {
		#A.ord = order(m$G>nash.G,m$sym.a,m$G,-(m$C-m$G),decreasing=TRUE)
		A.ord = order(!ignore.a,m$G>nash.G,m$sym.a,-(m$C-m$G),m$G,decreasing=TRUE)

	} else {
		A.ord = order(!ignore.a,m$G>nash.G,m$G,-(m$C-m$G),decreasing=TRUE)
	}
	m$G[A.ord]
	
	
	a.ind = 1
	a = A.ord[a.ind]
		
	while(a.ind <= length(A.ord)) {
		ret = get.Ue.L.a(m,a,lp=lp)
		if (is.null(ret$LY.mat)) {
			m$L[a] = Inf
			a.ind = a.ind +1
			a = A.ord[a.ind]
			next()
		} else {
			break()
		}
	}
			
	lp = ret$lp
	
	if (is.na(m$L[a])) {
		m$L[a] = ret$LY.mat[1,"L"]
	}	
	
	collab = c("L","Ue","a")
	Ue.opt = cbind(ret$LY.mat[,1],ret$LY.mat[,2],a)
	colnames(Ue.opt) = collab
	#rownames(Ue.opt)[1] = names(m$G)[a]
		
	Y.min = -Inf
	# Add the best Nash equilibrium to Ue.opt
	if (NROW(m$nash)>0 & Ue.opt[1,1]>0) {
		nash.G = max(m$G[m$nash])
		a = which(m$G == nash.G & get.is.nash(m))[1]
		LY.mat = cbind(0,nash.G,a)
		colnames(LY.mat) = collab
		
		Ue.opt = LY.get.upper.envelope(Xo=Ue.opt[,1],Yo=Ue.opt[,2],ao=Ue.opt[,3],
		                               Xn = LY.mat[,1], Yn = LY.mat[,2], an = a, collab=collab)
		                               
		if(sum(!is.finite(Ue.opt[,1]))>0) {
			warning("NA's in Ue.opt (adding Nash) This should not be. Debug!")
			stop()
		}
		                               		                               
   	Y.min   = nash.G                         
	}

	plot.main = paste("Ue (#A ", max.a.check, " / " ,m$nA, ")", sep="")
  if (do.plot) {		                             
	  plot(Ue.opt[,1],Ue.opt[,2],type="o",col="black",main=plot.main)
  }
  
	a.ind.start = a.ind +1
	recalc.opt.approx = TRUE
	num.L.calc = 2			
	for (a.ind in a.ind.start:m$nA) {
		#a.ind = a.ind +1
		if (a.ind>max.a.check) {break();}

		#flush.console()
		a = A.ord[a.ind]
		
		#if (a==27)
		#  restore.point("loc27")
    #restore.point("loc27")		  
    
    
		if (recalc.opt.approx) {
			if (sum(!is.na(Ue.opt[,1]))>1) {
				opt.approx = approxfun(Ue.opt[,1], y = Ue.opt[,2], method="linear",
	                      yleft=-Inf, yright=max(Ue.opt[,2]), rule = 1, f = 0, ties = "ordered")
			} else {
				opt.approx = approxfun(c(Ue.opt[1,1],Ue.opt[1,1]+1), y = c(Ue.opt[1,2],Ue.opt[1,2]), method="linear",
	                      yleft=-Inf, yright=max(Ue.opt[,2]), rule = 1, f = 0, ties = "ordered")
	  	}
  	}
		recalc.opt.approx = FALSE	
		# Check whether we can neglect the profile automatically
		if (!is.na(m$L[a])) {
			if (opt.approx(m$L[a]) >= m$G[a])	next();
		} else {
			if (opt.approx(m$C[a]-m$G[a]) >= m$G[a])	next();
		}
		num.L.calc = num.L.calc +1
		if (do.plot) {
  		try(points(m$C[a]-m$G[a],m$G[a], col="yellow"))
	  	try(legend("topleft",legend=paste("a:", a.ind,"L:",num.L.calc),bg="lightgrey"))
    }
		# If we cannot neglect it, let us do the whole calculation
		# SPEED UP THE ALGORITHM AT THIS POINT LATER!
		
		ret = get.Ue.L.a(m,a,lp=lp,Y.min = Y.min,LY.opt = Ue.opt, opt.approx = opt.approx)
		
		#plot(ret$LY.mat[,"L"],ret$LY.mat[,2],type="l")
		#ret1 = get.Ue.L.a(m,a,lp=lp)
		#lines(ret1$LY.mat[,"L"],ret1$LY.mat[,2],type="l",col="green")
		#lines(Ue.opt[,"L"],Ue.opt[,2],type="l",col="red",lty=2)
		
		
		recalc.opt.approx = ret$maybe.better.than.opt
		if (is.null(recalc.opt.approx)) {
			recalc.opt.approx = TRUE
		}
				
    if (is.null(ret$LY.mat)) {
			warning(paste("No L found for action profile ", a, " = ", m$lab.a[a]))
			print(paste("No L found for action profile ", a, " = ", m$lab.a[a]))
			m$L[a] = Inf
			next()
		}
    
    if (!is.null(m$L[a])) {
			m$L[a] = ret$LY.mat[1,"L"]
		}
		if (do.plot)
  		try(legend("topleft",legend=paste("a:", a.ind,"L:",num.L.calc),bg="lightgrey"))
		
		if (recalc.opt.approx) {
			Ue.opt = LY.get.upper.envelope(Xo=Ue.opt[,1],Yo=Ue.opt[,2],ao=Ue.opt[,3],
			                               Xn = ret$LY.mat[,1], Yn = ret$LY.mat[,2], an = a, collab=collab,
			                               plot.main = plot.main)
					            
	                          
			if(sum(!is.finite(Ue.opt[,1]))>0) {
				warning("NA's in Ue.opt This should not be. Debug!")
				stop()
			}
  	}
    #plot.Ue.opt(Ue.opt)
 		#a.ind = a.ind +1    	                               
	}
	
	if (sum(!is.na(Ue.opt[,1]))>1) {
		opt.approx = approxfun(Ue.opt[,1], y = Ue.opt[,2], method="linear",
                    yleft=-Inf, yright=max(Ue.opt[,2]), rule = 1, f = 0, ties = "ordered")
	} else {
		opt.approx = approxfun(c(Ue.opt[1,1],Ue.opt[1,1]+1), y = c(Ue.opt[1,2],Ue.opt[1,2]), method="linear",
                    yleft=-Inf, yright=max(Ue.opt[,2]), rule = 1, f = 0, ties = "ordered")
	}
  
	if (do.plot) {
	  plot.Ue.opt(Ue.opt,main=plot.main)
  }
	  
	rownames(Ue.opt) = rownames(m$g)[Ue.opt[,"a"]]
	Ue.opt = cbind(Ue.opt,m$G[Ue.opt[,"a"]]-Ue.opt[,"Ue"])
	colnames(Ue.opt)[NCOL(Ue.opt)] = "B"
	
	return(ret=list(m=m, opt.approx = opt.approx, Ue.opt = Ue.opt))
}

# Internal function
is.line.below.above.opt = function(L0,Y0,L1,Y1,LY.opt, opt.approx, below = TRUE) {
	
	#restore.point("is.line.below.above.opt")

	
	if (below) {
		# Check whether a point of LY.opt is below the line segment
		L.opt = which(LY.opt[,"L"] >= L0 & LY.opt[,"L"] <= L1 & !is.na(LY.opt[,"a"]))
		if (length(L.opt) > 0) {
			Y.line = 	( Y0 + ((Y1-Y0) / (L1-L0)) * (LY.opt[L.opt,"L"]-L0))
			if (sum(Y.line > LY.opt[L.opt,2]) > 0) {
				return(FALSE)
			}
		}
		# Check whether Y0 or Y1 are above LY.opt
		if (opt.approx(L0) < Y0 | opt.approx(L1) < Y1) {
			return(FALSE)
		}
		return(TRUE)
	} else {
		# Check whether a point of LY.opt is above the line segment
		L.opt = which(LY.opt[,"L"] >= L0 & LY.opt[,"L"] <= L1 & !is.na(LY.opt[,"a"]))
		if (length(L.opt) > 0) {
			Y.line = 	( Y0 + ((Y1-Y0) / (L1-L0)) * (LY.opt[L.opt,"L"]-L0))
			if (sum(Y.line < LY.opt[L.opt,2]) > 0) {
				return(FALSE)
			}
		}
		# Check whether Y0 or Y1 are below LY.opt
		if (opt.approx(L0) > Y0 | opt.approx(L1) > Y1) {
			return(FALSE)
		}
		return(TRUE)
	}		
}

# Internal function		
proceed.LY.mat = function(m,a=NULL,lp,LY.mat,LY.bound.mat,state="e",max.calc = 100,
                          row.LY = NULL, reduce = TRUE,tol= REP.CNT$TOL,Y.min = -Inf, Y.max = Inf,
                          LY.opt = NULL, opt.approx = NULL,use.LY.opt = !is.null(LY.opt)) {	

	
	#restore.point("proceed.LY.mat")
	
	
	if (state !="e") {
		i = state
		lambda.i = m$lambda[i]
	}
	
	# For the warning messages
	warn.lab.a = a
	G.a = m$G[a]


	
# 	
# 	plot(LY.mat,xlim=range(LY.mat[,1],LY.bound.mat[,1],na.rm=TRUE),
# 	            ylim=range(LY.mat[,2],LY.bound.mat[,2],na.rm=TRUE))
# 	lines(LY.mat)	            
# 	points(LY.bound.mat, col="blue")
# 	            	            
	STOP.calc = 10000
	nr = m$ny+1
	i.calc = 0
	if (is.null(row.LY)) {
		row.LY = max(which(!(is.na(LY.mat[,"L"]))))
	}
	row.low.LY = max(which(!(is.na(LY.bound.mat[,"L"]))))
	
	
	maybe.better.than.opt = !use.LY.opt
	
	while(sum(!is.na(LY.bound.mat[,1]) & LY.bound.mat[,"proceed"])>0) {
		#Check whether LYmat or LY.bound.mat have to be enlarged
		if (row.LY >= NROW(LY.mat)) {
			LY.mat = rbind(LY.mat,matrix(NA,nr,4))
		}
		if (row.low.LY >= NROW(LY.bound.mat)-1) {
			LY.bound.mat = rbind(LY.bound.mat,matrix(NA,nr,6))
		}
		
		# Get point in lowLY with highest gap
		rows = LY.bound.mat[,"proceed"] & (!is.na(LY.bound.mat[,1]))
		max.gap = max(abs(LY.bound.mat[rows,"gap"]),na.rm=TRUE)
		
		low.i = which(abs(LY.bound.mat[,"gap"])==max.gap & LY.bound.mat[,"proceed"])[1]
		L = LY.bound.mat[low.i,"L"]
		
		is.below = c(FALSE,FALSE)
		# Check whether the parts of the function to the left and right of point
		# low.i can be better than the optimal function calculated so far
		if (use.LY.opt) {
			below = (state == "e")
			i0 = LY.bound.mat[low.i,"i0"];i1 = LY.bound.mat[low.i,"i0"];
			is.below[1] = is.line.below.above.opt(L0=LY.mat[i0,"L"],Y0=LY.mat[i0,"Y"],
			                  L1=LY.bound.mat[low.i,"L"],Y1=LY.bound.mat[low.i,"Y"],LY.opt, opt.approx, below=below)
			is.below[2] = is.line.below.above.opt(L0=LY.bound.mat[low.i,"L"],Y0=LY.bound.mat[low.i,"Y"],
																			L1=LY.mat[i1,"L"],Y1=LY.mat[i1,"Y"],LY.opt, opt.approx,below=below)
		}
		
		# We cannot improve with this point, we can neglect it without calculation
		if (sum(is.below)==2) {
			# Remove the actual point from LY.bound.mat
			LY.bound.mat[low.i,]=NA
			next;
		} else {
			# Very rough indicator that the new function may improve on the envelope
			# At least the first upper envelope (with 3 points) could improve
			maybe.better.than.opt = TRUE
		}
		# Calculate U(L|a) for this point and add to LY.mat
		if (state=="e") {
			# Fix the resistance to the given value
			#glp_set_col_bnds(lp,1,GLP_FX,L,L)
	
			change.lp.get.B(lp=lp,m=m,a=a,L=L, only.L.change = TRUE)
			ret = solve.glpk.lp(lp,call.lab=paste("proceed.LY.mat: get.B a=",m$lab.a[a]," L = ", L,sep=""))
			Y = G.a-ret$objval
			slope = -lp.get.col.shadow.price(lp,1)
			
		# Calculate vi(L|a) for this point and add to Owi.mat
		} else {
			change.lp.get.pi.i(lp=lp,m=m,a=a,i=i,L=L)
			ret = solve.glpk.lp(lp,call.lab=paste("proceed.LY.mat: get.pi.i a=",m$lab.a[a]," L = ", L,sep=""))
			Y = m$g[a,state]+ret$objval  +lambda.i*L
			slope = lp.get.col.shadow.price(lp,1) + lambda.i
  	}
		row.LY = row.LY+1
		LY.mat[row.LY,] = c(L,Y,slope,TRUE)

#  		plot(LY.mat,xlim=range(LY.mat[,1],LY.bound.mat[,1],na.rm=TRUE),
#  	            ylim=range(LY.mat[,2],LY.bound.mat[,2],na.rm=TRUE))
#  		lines(LY.mat)	            
#  		points(LY.bound.mat, col="blue")
# 		for (i in 1:NROW(LY.mat)) {
# 			if (is.na(LY.mat[i,"L"])) {break()}
# 				draw.line.with.point.and.slope(LY.mat[i,"L"],LY.mat[i,"Y"],LY.mat[i,"slope"],
# 			                               lty=2,col="grey")
# 		}
# 		points(LY.bound.mat[low.i,"L"],LY.bound.mat[low.i,"Y"],col="yellow")
# 		points(LY.mat[row.LY,"L"],LY.mat[row.LY,"Y"],col="red")
		
		if (LY.mat[row.LY,1]==LY.mat[row.LY-1,1]) {
			warning("Repetitions in LY.mat")
			stop()
		}
		
		# We must check now whether we are at a kink. There is no correct slope in a kink
		# Add new lowLY only if it is no kink
		if (abs(Y-LY.bound.mat[low.i,"Y"]) > tol) {
			LY.mat[row.LY,"kink"] = FALSE
			# Create for each of the two neighbour segments of L 
			# a new point in LY.bound.mat

			############################################################################################			
			# Point for the segment to the left
			#############################################################################################
			
			if (!is.below[1]) {
				pi0=LY.bound.mat[low.i,"i0"];
				i0=pi0;i1=row.LY
				point = get.lb.L.Y(i0,i1, LY.mat,tol)
				
			#	points(point[1],point[2],col="orange")
				if (!is.null(point)) {
					if (point["X"] > LY.mat[i0,1] & point["X"] < LY.mat[i1,1]) {
						row.low.LY = row.low.LY + 1
						LY.bound.mat[row.low.LY,]=c(point,TRUE,i0,i1)
						if (LY.mat[i1,"Y"] < Y.min | LY.mat[i1,"Y"] > Y.max) {
							LY.bound.mat[row.low.LY,"proceed"]=FALSE
						}
					} else {
  					warning.proceed.LY.mat.point.outside(m,a,LY.mat,point,i0,i1)
					}
				}
				#LY.bound.mat[-1,] = NA
			}			
			############################################################
			# Point for the segment to the right
			############################################################
			if (!is.below[2]) {
				pi1=LY.bound.mat[low.i,"i1"];
				i0=row.LY;i1=pi1;
				point = get.lb.L.Y(i0,i1, LY.mat,tol)
			#	points(point[1],point[2],col="orange")
	
				if (!is.null(point)) {
					if (point["X"] > LY.mat[i0,1] & point["X"] < LY.mat[i1,1]) {
						row.low.LY = row.low.LY + 1
						LY.bound.mat[row.low.LY,]=c(point,TRUE,i0,i1)
						if (LY.mat[i1,"Y"] < Y.min | LY.mat[i1,"Y"] > Y.max) {
							LY.bound.mat[row.low.LY,"proceed"]=FALSE
						}
					} else {
  					warning.proceed.LY.mat.point.outside(m,a,LY.mat,point,i0,i1)
					}
				}
			}
		}
		# Remove the actual point from LY.bound.mat
		LY.bound.mat[low.i,]=NA
		#points(LY.bound.mat, col="orange")

						
		i.calc = i.calc +1
		if (i.calc >= STOP.calc) {
			stop("proceed.LY.mat hit STOP.calc please debug")	
		}
		
		
		if (i.calc >= max.calc) {
			warning("proceed.LY.mat max calc hit")
			break()
		}
		
		#if (i.calc == 6) break;

	} # End while
	
	if (reduce) {
		LY.mat = LY.mat[!is.na(LY.mat[,1]),,drop = FALSE]
		LY.mat = LY.mat[order(LY.mat[,"L"]),,drop = FALSE]
	}
	
	if (sum(!is.na(LY.bound.mat))==0) {
		LY.bound.mat = NULL
		max.gap = 0
		#LY.mat = LY.mat[LY.mat[,"kink"]==1,]
	} else {
		if (reduce) {
			LY.bound.mat = LY.bound.mat[!is.na(LY.bound.mat[,1]),,drop = FALSE]
			LY.bound.mat = LY.bound.mat[order(LY.bound.mat[,"L"]),,drop = FALSE]
		}
		max.gap = max(LY.bound.mat[,"gap"])
	}		
	
	return(list(LY.mat=LY.mat, LY.bound.mat = LY.bound.mat, max.gap=max.gap, lp=lp, maybe.better.than.opt=maybe.better.than.opt))
}

# Internal function
warning.proceed.LY.mat.point.outside = function(m,a,LY.mat,point,i0,i1) {
	str = paste("proceed.LY.mat a=",m$lab.a[a], " (Left segment) A point was outside the segment. Point ignored",
  					  "L[i0]=", LY.mat[i0,1], ",L.cut=", point["X"], ", L[i1]=", LY.mat[i1,1]) 
	#warning(str)
	#print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	print(paste(str))
}

# Get liquidity requirement of action profile a
get.L = function(m,a,lp=m$lp) {
	
	#assign.list(VAR.STORE[["get.L"]]$var)

	# Assume that we can implement every action profile without money burning
	make.ac.for.lp(m,a=a)
	change.lp.get.L(lp,m,a)
	ret = solve.glpk.lp(lp,call.lab=paste("get.L a=",m$lab.a[a],sep=""))
	if (ret$status != 0) {
		return(list(L=Inf,lp=lp))
	}
	return(list(L=ret$objval,lp=lp))
}	

# Get liquidity requirement and the amount of money burning neccessary at L(a)
# also get the minimal required expected amount of money burning to implement a
# and the corresponding liquidity reqiremen L.max(a)
get.L.and.B = function(m,a=NULL,lp=m$lp,sym.a=NULL,sym.y=NULL) {
	
	#restore.point("get.L.and.B")

	make.ac.for.lp(m,a=a)
	nr = m$lp.info$nr
	change.lp.get.L.of.B(lp=lp,m=m,a=a,B=0)
	ret = solve.glpk.lp(lp,retry.with.standard.basis = TRUE, warning.if.no.solution = FALSE, should.always.solve=FALSE,
	                    call.lab=paste("get.L.and.B get.L(B=0) a=",m$lab.a[a],sep=""))
		
	if (ret$status == 0) {
		B.min = 0
		B.min.slope = 1/ret$shadow.price[1]
		L.max = ret$objval
	# Zero Money burning cannot be implemented
	} else {
		change.lp.get.min.B(lp=lp,m=m,a=a)
		glp_std_basis(lp);
		ret = solve.glpk.lp(lp,call.lab=paste("get.L.and.B get.min.B a=",m$lab.a[a],sep=""))
		# Nothing can be implemented
		if (ret$status != 0) {
				return(list(L=Inf,lp=lp))
		}
		B.min = ret$objval

		# Get L from B
		change.lp.get.L.of.B(lp=lp,m=m,a=a,B=B.min)
		change.lp.get.L.of.B(lp=lp,m=m,a=a,B=B.min*1.1)
		ret = solve.glpk.lp(lp,
		      call.lab=paste("get.L.and.B: get.L(B.min) a=",m$lab.a[a]," B.min = ",B.min,sep=""))
		L.max = ret$objval
		B.min.slope = 1/ret$shadow.price[1]
	}
	
	change.lp.get.L(lp=lp,m=m,a=a)
	ret = solve.glpk.lp(lp,call.lab=paste("get.L.and.B: get.L a=",m$lab.a[a],sep=""))
	if (ret$status != 0) {
	  return(list(L=NA,lp=lp))
	}
	L.min = ret$objval
	
	if (L.min < L.max) {
		change.lp.get.B(lp=lp,m=m,a=a,L=L.min*(1+(10^(-7))))
		ret = solve.glpk.lp(lp,
		      call.lab=paste("get.L.and.B get.B(L.min) a=",m$lab.a[a]," L.min = ",L.min,sep=""))
		      
		if (ret$status != 0) {
  		L.min = L.min*(1+(10^(-7)))
  		
  		change.lp.get.B(lp=lp,m=m,a=a,L=L.min)
		  ret = solve.glpk.lp(lp,
		      call.lab=paste("get.L.and.B get.B(L.min) a=",m$lab.a[a]," L.min = ",L.min,sep=""))
  		if (ret$status == 0) {
    		str = paste("get.L.and.B get.B(L.min) a=",m$lab.a[a]," L.min = ",L.min, 
    		      " only solved after L.min was multiplied by 1+(10^(-7))", sep="")
        print(str)
        warning(str)
      } else {		      
    		str = paste("get.L.and.B get.B(L.min) a=",m$lab.a[a]," L.min = ",L.min, 
    		      " could not be solved even after L.min was multiplied by 1+(10^(-7))", sep="")
        print(str)
        warning(str)
        return(list(L=NA,lp=lp))
      }
  	}
		B.max = ret$objval
		B.max.slope = lp.get.col.shadow.price(lp,1)
	}
	if (L.min >= L.max) {
		L.max = L.min
		B.max = B.min
		B.max.slope = NA
		B.min.slope = NA
	}
	return(list(L=L.min, L.max = L.max, B.max = B.max, B.min=B.min,
	            B.max.slope= B.max.slope, B.min.slope = B.min.slope,lp=lp))
}	


# Implement actions without money burning
get.Ub = function(m,a,L=Inf) {
  
	#restore.point("get.Ub")

	# Assume that we can implement every action profile without money burning
	make.ac.for.lp(m,a=a)
	change.lp.max.exp.B(m$lp,m,a,L)
	ret = solve.glpk.lp(lp,call.lab=paste("get.Ub a=",m$lab.a[a],sep=""))
	if (ret$status != 0) {
		return(list(Ub=NA))
	}
	return(list(Ub=ret$objval))
}

#######################################################################################################
# Several functions that
# set the linear program to a particular task
#######################################################################################################

change.lp.max.exp.B = function(lp,m,a,L) {
  # Assign phi.mat for a
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}

  
	# Fix the resistance to the given value
	glp_set_col_bnds(lp,1,GLP_FX,L,L)

	lab.a = make.lab.a.ai.prob(m=m,a=a)
	lp.set.name(lp, paste("get.max.exp.B(L=", round(L,5)," | a=",lab.a,")", sep="" ))
	
	# Unfix B
	glp_set_row_bnds(lp,1,GLP_FR,0,0)

	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}
	

	# Simply the opposite sign from the problem of minimizing expected B
	obj.coef = -make.obj.coef.for.B(m=m,a=a)
	lp.set.obj.coef(lp,obj.coef)
}


change.lp.get.pi.i = function(lp,m,a,i,L,only.L.change = FALSE) {
	if (is.null(a)) {
		stop("change.lp.get.pi.i not yet implemented for mixed strategies")
	}

	# Assign phi.mat for a
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}
	
	lp.set.name(lp, paste("get.pi.i",i,",(L=", round(L,5)," | a=",a,")", sep="" ))
	
	# Fix the resistance to the given value
	glp_set_col_bnds(lp,1,GLP_FX,L,L)
		
	if (only.L.change) return();

	# Expected payments of player i. They get negative weight, since we have a minimzation problem
	obj.coef = rep(0,m$n*m$ny+1)
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}
	obj.coef[(2+(i-1)*m$ny):(1+i*m$ny)] = - m$phi.mat[,a]
	lp.set.obj.coef(lp,obj.coef)
	
	# Unfix B
	glp_set_row_bnds(lp,1,GLP_FR,0,0)
}

change.lp.get.B = function(lp,m,a,L,only.L.change = FALSE) {
  # Assign phi.mat for a
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}

  
	# Fix the resistance to the given value
	glp_set_col_bnds(lp,1,GLP_FX,L,L)
	
	if (only.L.change) return();

	lab.a = make.lab.a.ai.prob(m=m,a=a)
	lp.set.name(lp, paste("get.B(L=", round(L,5)," | a=",lab.a,")", sep="" ))
	
	# Unfix B
	glp_set_row_bnds(lp,1,GLP_FR,0,0)
	
	obj.coef = make.obj.coef.for.B(m=m,a=a)
	lp.set.obj.coef(lp,obj.coef)
}

change.lp.get.L = function(lp,m,a) {
  
  # Assign phi.mat for a
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}
	
	lp.set.name(lp, paste("get.L(a=",a,")", sep="" ))
	
	# Allow the resistance to take any non-negative value
	glp_set_col_bnds(lp,1,GLP_LO,0,0)

	# Unfix B
	glp_set_row_bnds(lp,1,GLP_FR,0,0)

	# Label		
	lab.a = make.lab.a.ai.prob(m=m,a=a)
	lp.set.name(lp, paste("get.L(a=",lab.a,")", sep="" ))

	# Objective function
	obj.coef = make.obj.coef.for.L(m=m,a=a)
	lp.set.obj.coef(lp,obj.coef)
}	

change.lp.get.L.of.B = function(lp,m,a,B) {
	
	#assign.list(VAR.STORE[["change.lp.get.L.of.B"]]$var)

  # Assign phi.mat for a
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}
	
			
	# Allow the resistance to take any non-negative value
	glp_set_col_bnds(lp,1,GLP_LO,0,0)

	# Fix upper bound of B
	glp_set_row_bnds(lp,1,GLP_UP,0,B)
			
	lab.a = make.lab.a.ai.prob(m=m,a=a)
	lp.set.name(lp, paste("get.L(B=", round(B,5)," | a=",lab.a,")", sep="" ))

	obj.coef = make.obj.coef.for.L(m=m,a=a)
	lp.set.obj.coef(lp,obj.coef)
}

change.lp.get.min.B = function(lp,m,a=NULL) {

  # Assign phi.mat for a
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[a]]}
	
	# Allow the resistance to take any non-negative value
	glp_set_col_bnds(lp,1,GLP_LO,0,0)

	# Unfix B
	glp_set_row_bnds(lp,1,GLP_FR,0,0)
		
	
	
	lab.a = make.lab.a.ai.prob(m=m,a=a)
	lp.set.name(lp, paste("get.min.B(a=",lab.a,")", sep="" ))

	obj.coef = make.obj.coef.for.B(m=m,a=a)
	lp.set.obj.coef(lp,obj.coef)
}


########################################################################################################################
# Functions for initialization of LP for a given a
########################################################################################################################

make.lab.a.ai.prob= function(m,a=NULL,ai.prob=NULL,private=FALSE) {
	if (!is.null(a)) {
		lab = paste(a,m$lab.a[a])
	} else if (!is.null(ai.prob)) {
		if (private) {
			lab = "mixed priv"
		} else {
			lab = "mixed pub"
		}
	}
	lab
}
		

make.obj.coef.for.B = function(m,a=NULL,ai.prob=NULL,private=FALSE) {
	if (!is.null(a)) {
		obj.coef = c(0,rep(m$phi.mat[,a], times=m$n))
	} else if (!is.null(ai.prob)) {
		if (private) {
			p.y.a.coef = get.p.y.a.coef(m,ai.prob)
			obj.coef = c(0,rep(p.y.a.coef, times=m$n))
		} else {
			p.y.coef = get.p.y.coef(m,ai.prob)
			obj.coef = c(0,rep(p.y.coef, times=m$n))
		}
	}
	return(obj.coef)
}

make.obj.coef.for.L = function(m,a=NULL,ai.prob=NULL,private=FALSE) {
	if (!is.null(a)) {
			# Set objective to the resistance L
			obj.coef = c(1,rep(0, length=m$n*m$ny))
	} else if (!is.null(ai.prob)) {
		if (private) {
			obj.coef = c(1,rep(0, length=m$n*m$ny*m$nA))
		} else {
			obj.coef = c(1,rep(0, length=m$n*m$ny))
		}
	}
	return(obj.coef)
}
	
make.constraints.for.a = function(m,ak, old.con.struct = m$old.con.struct, L.weights = m$L.weights) {
	
	#assign.list(VAR.STORE[["make.constraints.for.a"]]$var)		
	
	if (is.null(L.weights)) {
		L.weights = rep(1/m$n,m$n)
	}
	
	n = m$n; ny = m$ny; a.dim = m$a.dim;
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[ak]]}
	phi.ak = m$phi.mat[,ak]
	g.ak = m$g[ak,]
	
	
	all.new = is.null(old.con.struct)

	if (all.new) {		
		ncon = n*ny + sum(m$a.dim) - m$n + m$ny +1
		con = matrix(0,nrow = ncon, ncol = 1+ny*n)
		rhs = rep(NA,ncon)
		dir = rep(NA,ncon)
		bounds = list(lower=rep(-Inf,NCOL(con)),upper=rep(Inf,NCOL(con)))
		
		str = NULL
		for (i in 1:m$n) {
			str = c(str,paste("p",i,"(",rownames(m$phi.mat),")",sep=""))
	  }
		colnames(con)=c("L",str)
		rownames(con)=1:NROW(con)
	} else {
		con = old.con.struct$con
		rhs = old.con.struct$rhs
		ncon = NROW(con)
	}
	
	# Payment IC for every player and every signal
	row = 0
	for (i in 1:n) {
		for (y in 1:ny) {
			row = row+1
			con.row = matrix(0,ny,n)	
			con.row[y,i] = 1
			con[row,] = c(-L.weights[i],as.vector(con.row))
		}
		if (all.new) {
			rownames(con)[(row-ny+1):row] = paste("PC p",i,"(",m$lab.y[1:ny],")",sep="")
		}
	}

	if (all.new) {
		rhs[1:row] = 0
		dir[1:row] = "<="
	}
	#con
	
	# Action IC for every player i and all possible action profiles
	rowstart = row +1
	for (i in 1:n) {
		replies.i = get.replies.ind(m,a.ind=ak,i=i)
		A.act = setdiff(replies.i,ak)
		for (a in A.act) {
			con.row = matrix(0,ny,n)
			con.row[,i] = m$phi.mat[,a]-m$phi.mat[,ak]
			row = row+1
			con[row,] = c(0,as.vector(con.row))
			rhs[row] = m$g[a,i]-m$g[ak,i]
		}
		if (all.new) {
			rownames(con)[(row-NROW(A.act)+1):row] = paste("AC i",i," ",m$lab.a[A.act],sep="")
		}
	}
	if (all.new) {
		dir[rowstart:row]=">="
	}
	#con

	if (all.new) {		
		# Budget Constraint for every signal y
		rowstart = row +1
		for (y in 1:ny) {
			con.row = matrix(0,ny,n)
			con.row[y,] = 1
			row = row+1
			con[row,] = c(0,as.vector(con.row))
		}
		rownames(con)[rowstart:row] = paste("BC ",m$lab.y[1:ny],sep="")
	
		rhs[rowstart:row]=0
		dir[rowstart:row] = ">="
		
		# Add the total money burning constraint
		# Even though with B=0, we could incorporate this constraint
		# via the budget constraints, using this constraint
		# allows to retrieve the slope of B(L|a) at B=0 via the shadow price
		# also we need to adapt this constraint in case B=0 is never feasible
		row = row+1
		con[row,] = c(0,rep(m$phi.mat[,ak], times=m$n))
		rownames(con)[row] = "Bmax"
		rhs[row] = 0
	  dir[row] = "free"
		ret = list(obj.coef = NULL, con = con, dir = dir, rhs = rhs, obj.max = FALSE, 
		           bounds = bounds, row.lab=rownames(con), col.lab=colnames(con))	                         
	} else {
		# B constraint
		con[NROW(con),] = c(0,rep(m$phi.mat[,ak], times=m$n))
		ret = list(obj.coef = NULL, con = con, dir = old.con.struct$dir, rhs = rhs, obj.max = FALSE, 
		           bounds = old.con.struct$bounds, row.lab=rownames(con), col.lab=colnames(con))	                         
	}
		
	return(ret)
}

make.lp.for.a = function(m,ak, old.lp, L.weights = m$L.weights) {
	
	#assign.list(VAR.STORE[["make.lp.for.a"]]$var)		
		
	if (is.null(L.weights)) {
		L.weights = rep(1/m$n,m$n)
	}
	
	n = m$n; ny = m$ny; a.dim = m$a.dim;
	if (!is.null(m$phi.mat.li)) {m$phi.mat = m$phi.mat.li[[ak]]}
	phi.ak = m$phi.mat[,ak]
	g.ak = m$g[ak,]
	
	npc = n*ny; nac =  sum(m$a.dim)- m$n; nbc = m$ny;
	start.pc = 1; start.ac = npc+1; start.bc = npc+nac+1;
	
	nr = npc + nac +nbc;
	nc = 1+ny*n
	
	if (is.null(old.lp)) {
		lp = create.lp.prob(nr,nc)
		
		str = NULL
		for (i in 1:m$n) {
			str = c(str,paste("p",i,"(",rownames(m$phi.mat),")",sep=""))
	  }
		collab = c("L",str)
		rowlab = rep("",nr)

		new.lp = TRUE
	} else {
		new.lp = FALSE
		lp = old.lp
	}
	
	# Set payment IC for every player and every signal
	if (new.lp) {	
		row = start.pc-1
		for (i in 1:n) {
			for (y in 1:ny) {
				row = row+1
				set.row(lp, row, c(-L.weights[i],1), indices = c(1,(i-1)*ny+y+1))
			}
			rowlab[(row-ny+1):row] = paste("PC p",i,"(",m$lab.y[1:ny],")",sep="")
		}
		set.rhs(lp,rep(0,npc),start.pc:(start.pc+npc-1))
		set.constr.type(lp, rep("<=",npc), start.pc:(start.pc+npc-1))
	}

		
	# Action IC for every player i and all possible action profiles
	row = start.ac-1
	for (i in 1:n) {
		replies.i = get.replies.ind(m,a.ind=ak,i=i)
		A.act = setdiff(replies.i,ak)
		for (a in A.act) {
			con.row = matrix(0,ny,n)
			con.row[,i] = m$phi.mat[,a]-m$phi.mat[,ak]
			row = row+1
			set.row(lp,row,c(0,as.vector(con.row)))
			set.rhs(lp,m$g[a,i]-m$g[ak,i],row)
		}
		#rowlab[(row-NROW(A.act)+1):row] = paste("AC i",i,sep="")
		rowlab[(row-NROW(A.act)+1):row] = paste("AC i",i," ",m$lab.a[A.act],sep="")
		
	}
	if (new.lp) {
		set.constr.type(lp, rep(">=",nac),start.ac:(start.ac+nac-1))
	}
		
	# Budget Constraint for every signal y
	if (new.lp) {		
		row = start.bc-1
		ind = ((1:n)-1)*ny
		for (y in 1:ny) {
			row = row+1
			set.row(lp,row,rep(1,n),ind+y+1)
		}
		set.rhs(lp,rep(0,nbc),start.bc:(start.bc+nbc-1))
		set.constr.type(lp, rep(">=",nbc), start.bc:(start.bc+nbc-1))
		rowlab[start.bc:(start.bc+nbc-1)] = paste("BC ",m$lab.y[1:ny],sep="")
	}
	
	if (new.lp) {
		set.bounds(lp,c(0,rep(-Inf,nc-1)),rep(Inf,nc))
		dimnames(lp) <- list(rowlab,collab)
	} else {
		dimnames(lp)[[1]][start.ac:(start.ac+nac-1)] <- rowlab[start.ac:(start.ac+nac-1)]
	}
	return(lp)
}



#######################################################################################################################
# Running LP-OPS
#######################################################################################################################

# Run the big optimization problem
run.LP.OPS = function(m,ae=1, ai=c(4,4),delta=0.5, return.details = TRUE,ae.sym.pay.signals=NULL, target = "Ue-V") {
		
  #restore.point("run.LP.OPS")		
	
	n = m$n; ny = m$ny; a.dim = m$a.dim;

	if ( (!m$sym.a[ae]) &  !is.null(ae.sym.pay.signals)) {
		warning(paste("run.LP.OPS: ae.sym.pay.signals was set to null since ae is not symmetric"))
		ae.sym.pay.signals = NULL
	}
	
	phi.mat.k = list()
	if (!is.null(m$phi.mat.li)) {
		phi.mat.k$e = m$phi.mat.li[[ae]]
		for (i in 1:n) {
			phi.mat.k[[i+1]] = m$phi.mat.li[[ai[i]]]
		}
	} else {
		phi.mat.k$e = m$phi.mat
		for (i in 1:n) {
			phi.mat.k[[i+1]] = m$phi.mat
		}
	}
	
	if (!is.numeric(ae)) {
		ae = get.a.by.name(m,ae)
	}
	if (!is.numeric(ai)) {
		ai = get.a.by.name(m,ai)
	}

	as = c(ae,ai)
	
	

		
  
	ncon = (n*ny + sum(a.dim) - n + ny)*(n+1)
	con = matrix(0,nrow = ncon, ncol = ny*n * (n+1))
	rhs = rep(NA,ncon)
	dir = rep(NA,ncon)
	
	str = NULL
	for (k in 1:(n+1)) {
		for (i in 1:n) {
			str = c(str,paste("p",c("e",1:n)[k],".",i,"(",rownames(phi.mat.k[[k]]),")",sep=""))
  	}
	}
	colnames(con)=str
	rownames(con)=1:NROW(con)

	# Goal maximize Ue-V.
	# Since Ue = G(ae)-Sum[pe] and V=sum(ci(ai))+sum(pi[ai])
	# we have to maximize
	# -(Sum(pe)-Sum(pi.i[ai]))
	
	obj.coef = rep(0,NCOL(con))
	if (target == "Ue-V" | target == "Ue") {
		obj.coef = c(-rep(phi.mat.k[["e"]][,ae], times=n),rep(0,times=ny*n*n))  
	} 
	if (target == "V" | target == "Ue-V") {
		for (i in 1:n) {
			start = (ny*n*(i)+ny*(i-1))+1
			obj.coef[start:(start+ny-1)] = phi.mat.k[[i+1]][,ai[i]]
	  }
	}
	# minimize vi
	if (is.numeric(target)) {
			i = target
			start = (ny*n*(i)+ny*(i-1))+1
			obj.coef[start:(start+ny-1)] = phi.mat.k[[i+1]][,ai[i]]
	}
	names(obj.coef) = str
	obj.coef		
	
	k = 1:(n+1)
	start.k = ny*n*(k-1) +1
	end.k   = start.k-1 +ny*n

	# Payment IC for every player, every signal y and every k
	#pik_hat(y)+delta*pie-delta*p11=delta(gi(ae)-g1(a1))
	row = 0
	row.start = row
	
	# PC for all k,i,y:
	# (1-delta)*pk.i(y)+delta*pe.i.exp-delta*pi.i.exp=delta*(g[ae,i]-g[ai,i])
	for (k in 1:(n+1)) {
		for (i in 1:n) {
			row.start.i = row
			for (y in 1:ny) {
				row = row+1
				# Equilibrium payments that influence ue
				con[row,] = 0
				con.row.ae = matrix(0,ny,n)	
				con.row.ae[,i] = delta * phi.mat.k[[1]][,ae]
				con[row,start.k[1]:end.k[1]] = as.vector(con.row.ae)
				
				# Punishment payments that influence vi
				con.row.ai = matrix(0,ny,n)	
				con.row.ai[,i] = -delta *  phi.mat.k[[i+1]][,ai[i]]
				con[row,start.k[i+1]:end.k[i+1]] =  con[row,start.k[i+1]:end.k[i+1]] + as.vector(con.row.ai)
				
				# Actual payments
				con.row = matrix(0,ny,n)	
				con.row[y,i] = 1
				con[row,start.k[k]:end.k[k]] = as.vector(con.row) + con[row,start.k[k]:end.k[k]] 
				rownames(con)[row] = paste("PC p",c("e",1:n)[k],".",i,"(",m$lab.y[y],")",sep="")
			}
			rhs[(row.start.i+1):row] = delta*(m$g[ae,i]-m$g[ai[i],i])
		}
	}
	dir[(row.start+1):row] = "<="
	con[row.start:row,]

	
	
	# Action IC for every k, every player i and every action of him
	rowstart = row +1
	for (k in 1:(n+1)) {
		m$phi.mat = phi.mat.k[[k]]
		for (i in 1:n) {
			replies.i = get.replies.ind(m,a.ind=as[k],i=i)
			for (a in setdiff(replies.i,as[k])) {
				con.row = matrix(0,ny,n)
				con.row[,i] = m$phi.mat[,a]-m$phi.mat[,as[k]]
				row = row+1
				con[row,start.k[k]:end.k[k]] = as.vector(con.row)
				rhs[row] = m$g[a,i]-m$g[as[k],i]
				#rownames(con)[row] = paste("ActIC k", k-1," i",i," ",m$lab.a[a],sep="")
				rownames(con)[row] = paste("AC a",c("e",1:n)[k],".",i," to ",m$lab.a[a],sep="")
			}
		}
	}
	dir[rowstart:row]=">="
	con
		
	# Budget Constraint in every state for every signal y
	rowstart = row +1
	for (k in 1:(n+1)) {
		for (y in 1:ny) {
			con.row = matrix(0,ny,n)
			con.row[y,] = 1
			row = row+1
			con[row,start.k[k]:end.k[k]] = as.vector(con.row)
			#rownames(con)[row] = paste("BC k",k-1,", y",y,sep="")
			rownames(con)[row] = paste("BC p.",c("e",1:n)[k],"(",m$lab.y[y],")",sep="")
		}
	}
	rhs[rowstart:row]=0
	dir[rowstart:row]=">="
#	con
#	cbind(con,dir,rhs)

	
	if (!is.null(ae.sym.pay.signals) & length(ae.sym.pay.signals)>0) {
		for (ind.y in ae.sym.pay.signals) {
			if (is.numeric(ind.y)) {
				y = ind.y
			} else {
				y = which(m$lab.y==ind.y)
			}
			for (i in 2:m$n) {
				con.row = rep(0,NCOL(con))
				con.row[c(y,ny*(i-1)+y)]=c(1,-1)
				con = rbind(con,con.row)
				dir = c(dir,"==")
				rhs = c(rhs,0)
				rownames(con)[NROW(con)]=paste("Sym pe(",m$lab.y[y],")",sep="")
			}
		}
	}
	
# 	for (i in 1:2) {
# 		con.row = rep(0,NCOL(con))
# 		con.row[ny*(i-1)+1]=1
# 		con = rbind(con,con.row)
# 		dir = c(dir,"==")
# 		rhs = c(rhs,0)
# 		rownames(con)[NROW(con)]=paste("pe(yC)=0",sep="")
# 	}
# 	
	
		
	apply(con,1,sum)
	
	bounds <- list(lower = list(ind = c(1:NCOL(con)), val = rep(-Inf,NCOL(con))),
	               upper = list(ind = c(1:NCOL(con)), val = rep(Inf,NCOL(con))))
	               

	lp = make.glpk.lp(obj.coef, con, dir, rhs, obj.max = TRUE, 
		             bounds, row.lab=rownames(con), col.lab=colnames(con), prob.name="LP-OPS", prob = NULL, lp=NULL,
		             set.names = TRUE)
	              
	ret = solve.glpk.lp(lp)
	#glp_print_sol(lp, "lp_sol_run.txt")	
	#glp_print_prob(lp, "lp_prob_run.txt")	

	
	if (return.details) {
		res   = list()
		p     = list()
		p.exp = matrix(NA,n+1,n)
		rownames(p.exp) = m$lab.a[as]
		for (k in 1:(n+1)) {
			p[[k]] = matrix(ret$solution[(ny*n*(k-1)+1):(ny*n*k)],ny,n)
			rownames(p[[k]]) = m$lab.y
			for (i in 1:n) {
				p.exp[k,i] = sum(p[[k]][,i] * phi.mat.k[[k]][,as[k]])
			}
		}
		names(p)=paste("k=",c("e",1:n),sep="")
		res$delta = delta
		res$p     = p
		res$p.exp = p.exp
		
		res$ue = m$g[ae,] - res$p.exp[1,]
		res$v  = (1-delta)*(m$g[cbind(ai,1:n)] - res$p.exp[cbind((2:(n+1)),1:n)]) + delta*res$ue
		res$p.max[1:n] = (delta / (1-delta)) * (res$ue-res$v)
		
		res$V = sum(res$v)
		res$Ue = sum(res$ue)
		res$B = as.numeric(m$G[ae]-res$Ue)
		res$status = ret$status
	}
	
	if (return.details) {
		return(res)
	} else {
		return(ret$opt)
	}
}	

# Gets the expected payments p from a vector of all p.hat
get.p.from.p.hat = function(m,as,p.hat) {
	ret = numeric(n)
	if (is.null(m$phi.mat.li)) {
		for (i in 1:n) {
			ret[i] = sum(p.hat*m$phi.mat[,as])
		}
	} else {
		for (i in 1:n) {
			ret[i] = sum(p.hat*m$phi.mat.li[[as]][,as])
		}
	}		
	return(ret)
}
get.p.from.p.hat = Vectorize(get.p.from.p.hat, vectorize.args = "as", SIMPLIFY = TRUE,USE.NAMES = TRUE)


make.constraints = function(m,a=NULL,ai.prob=NULL,private=TRUE, old.con.struct = m$old.con.struct, 
full.mix=FALSE, sym.y = NULL, sym.a = NULL) {
	if (!is.null(a)) {
		return(make.constraints.for.a(m=m,a=a,old.con.struct=old.con.struct))
	} else if (!is.null(ai.prob)) {
		if (private) {
			return(make.constraints.for.s.priv(m=m,ai.prob=ai.prob, old.con.struct = old.con.struct, 
                                  full.mix=full.mix, sym.y = sym.y, sym.a = sym.a))
   	} else {
	   	return(make.constraints.for.s.pub(m=m,ai.prob=ai.prob, old.con.struct = old.con.struct, 
                                  full.mix=full.mix, sym.y = sym.y))
  	}
	}
	warning("make.constraints must specify a (pure equilibria) or ai.prob (mixed equilibria)")
	return(NULL)
}	











