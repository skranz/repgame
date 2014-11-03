###################################################################################################################
###  Search functions to access certain action profiles
###################################################################################################################


get.is.nash = function(m) {
  m$is.nash = (m$C == m$G)
}
  
### Indeces of action profiles given labels of action profiles
get.a.by.name = function(m,a) {
	match(a,m$lab.a)
}


### will deprecat since 
get.a.vec = function(m,a.ind) {
	if (length(a.ind)>1) {
		return(t(sapply(a.ind,get.a.vec,m=m)))
	}
 	ceiling( (((a.ind-1) %% m$modulus.ai)+1) / m$shift.ai)
}


get.symmetric.a = function(m,a,source.i=1,dest.i=2) {
	if (m$n > 2) {
		warning("get.symmetric.a not yet implemented for more than 2 players!")
	}
	ai = get.ai.mat(m,a)
	cols = c(2:1)
	ai = ai[,cols]
	get.a.from.ai.mat(m,ai)
}

#a = 1:m$nA
#get.symmetric.a(m,1:m$nA)
	
### Gives the index of a specific action profile
get.a.ind = function(m,a.vec) {
	sum((a.vec-1)*m$shift.ai)+1
}

### a.ind is a matrix with each row characterizing an action profile
### in the form (a.1,a.2,...,a.n)
### get.a.vec returns a vector of internal indices for every such action profile 
get.a.from.ai.mat = function(m,ai.mat) {
	a = rep(1,NROW(ai.mat))
	for (i in 1:NCOL(ai.mat)) {
		a = a + (ai.mat[,i]-1)*m$shift.ai[i]
	}
	a
}

### Inverse function of get.a.from.ai.mat
get.ai.mat = function(m,a.ind=1:m$nA) {
	ai.mat = matrix(0,NROW(a.ind),m$n)
	for (i in 1:m$n) {
		ai.mat[,i] = ceiling( (((a.ind-1) %% m$modulus.ai[i])+1) / m$shift.ai[i])
	}
	ai.mat
}


### Calculates only column i of get.ai.mat
get.ai = function(m,a.ind,i) {	
 	ceiling( (((a.ind-1) %% m$modulus.ai[i])+1) / m$shift.ai[i])
}

### Get matrix of indices of those action profiles that are best replies
### given the action profiles in a.ind
get.br.mat = function(m,a.ind=NULL) {
  nr = m$nA
  if (!is.null(a.ind)) nr = NROW(a.ind)
  mat = matrix(NA,nr,m$n)
  for (i in 1:m$n) 
    mat[,i] = get.best.reply.ind(m,a.ind,i)
  mat
}

### Gives the indices of possible replies of player i holding fixed the actions of the other players
get.replies.ind = function(m,a.ind=NULL,i,a.vec=NULL) {
	if (is.null(a.ind)) {
		a.ind = get.a.ind(m,a.vec)
		ai = a.vec[i]
  } else {
		ai = get.ai(m,a.ind,i)
	}
	return ((a.ind - (ai-1)*m$shift.ai[i]) + ((1:m$a.dim[i])-1)*m$shift.ai[i])
}	  

get.best.reply.ind= function(m,a.ind=NULL,i) {
	replies = get.replies.ind(m,a.ind,i)
	return(replies[which.max(m$g[replies,i])])
}	  

is.unilateral.deviation = function(m, a, a.org) {
	if (NROW(a)>1) {
		return(mapply(is.unilateral.deviation, a = a, MoreArgs=list(m=m,a.org=a.org)))
	}
	a.vec     = get.a.vec(m,a)
	a.org.vec = get.a.vec(m,a.org)
	if (sum(a.vec != a.org.vec) != 1) {
		return(0)
	} else {
		return(which(a.vec != a.org.vec))
	}
}	


#' Returns m$opt.mat expanded to a vector of discount factors	
get.mat.of.delta = function(m,delta=NULL,r=NULL,add.crit.delta=TRUE, before.delta.gap=10^(-14), pm.keep.sparse = is.null(delta)) {
	
	#restore.point("get.mat.of.delta")

	if (is.null(delta)) {
		delta = 1 / 1+r
	}
	
	### Add delta at kinks and directly before kinks
	### This allows sufficiently precise plots
	if (add.crit.delta & !pm.keep.sparse) {
		rows = m$opt.mat[,"opt"] == 1 & m$opt.mat[,"delta"] >= min(delta) & m$opt.mat[,"delta"] <= max(delta)
		delta.add = m$opt.mat[rows,"delta"]
		delta = c(delta,delta.add-before.delta.gap,delta.add)
		delta = unique(sort(delta))
		delta = delta[delta>=0 & delta <1]
	}	
	r = (1-delta) / delta

		
	### Dealing with games of perfect monitoring
	if (m$sol.type != "imp") {
  	mat = m$opt.mat
  	mat = mat[mat[,"opt"] == 1,,drop=FALSE]
  	if (!pm.keep.sparse) {
    	min.delta = min(mat[,"delta"])
  	  rows = delta >= min.delta
  	  delta = delta[rows]
  	  r = r[rows] 

    	ind = findInterval(delta,mat[,"delta"])
    	mat = mat[ind,,drop=FALSE]
    	mat[,"delta"] = delta; mat[,"r"]  = r;
    	return(mat)
	  } else {
    	mat = pm.make.plot.mat(mat)
	  }  
  	return(mat)
	}
	
	ret = get.L.of.delta(m,delta,r=r,also.row.ind = TRUE)
	L = ret$L
	ind = ret$ind

	mat = m$opt.mat[ind,]
	mat[,"L"] = L; mat[,"delta"] = delta; mat[,"r"]  = r;
	mat[,"Ue"] = m$Ue.fun(L); mat[,"V"] = m$V.fun(L);
	mat[,"UV"] = mat[,"Ue"]-mat[,"V"]
	
	return(mat)	
}

#' Returns m$opt.mat expanded to a vector of liquidities
get.mat.of.L = function(m,L=seq(0,max(m$opt.mat[,"L"]),length=10000), add.crit.L = TRUE, eps = 10^(-9)) {
	
	#restore.point("get.mat.of.L")

	if ("helper" %in% colnames(m$opt.mat)) {	
    rows = m$opt.mat[,"helper"] == 0
  } else {
    rows = 1:NROW(m$opt.mat)
  }
	
  mat = m$opt.mat[rows,,drop=FALSE]	
	if (add.crit.L) {
  	L = sort(unique(c(L,mat[,"L"],mat[mat[,"opt"]==0,"L"]-eps)))
	}
	min.L = min(mat[,"L"])
	L = L[L>=min.L]
	
	ind = findInterval(L,mat[,"L"])
	mat = mat[ind,,drop=FALSE]
	colnames(mat)=colnames(m$opt.mat)
	mat[,"L"] = L
	if (m$sol.type=="imp") {
	  mat[,"Ue"] = m$Ue.fun(L);
	  mat[,"V"] = m$V.fun(L);
	  mat[,"UV"] = mat[,"Ue"]-mat[,"V"]
	  i.max = m$n
	  if (m$symmetric) i.max = 1
	  for (i in 1:i.max)
	    mat[,paste("v",i,sep="")] = m$vi.fun[[i]](L);
  }
	return(mat)	
}

#' Gets the vector of the maximum liquidities that can be generated for given 
get.L.of.delta = function(m,delta=NULL,r=NULL, also.row.ind = FALSE) {
	
	#restore.point("get.L.of.delta")
		
	opt.mat = m$opt.mat
	if (is.null(delta)) {
		delta = 1 / 1+r
	}
	
	if ("helper" %in% colnames(m$opt.mat)) {	
  	rows = which(opt.mat[,"helper"] == 0 & opt.mat[,"opt"] == 1)
  } else {
    rows = which(opt.mat[,"opt"] == 1)
  }
	ind = findInterval(delta,opt.mat[rows,"delta"])
	ind = rows[ind]
	L = opt.mat[ind,"L"]
	
	exact = opt.mat[ind,"delta"] == delta | !opt.mat[ind,"opt.next"]
	### For points on a line, we use the formula
	### L(r) =	((L1y2-L2y1)/(y2-y1-r(L2-L1)))
	or = ind[!exact]
#	interc = ( (opt.mat[or,"L.next"]*opt.mat[or,"UV"]) -  (opt.mat[or,"L"]*opt.mat[or,"UV.next"]) ) /
						#(opt.mat[or,"L.next"]-opt.mat[or,"L"])
	
	L[!exact] = ( (opt.mat[or,"L"]*opt.mat[or,"UV.next"]) - (opt.mat[or,"L.next"]*opt.mat[or,"UV"]) ) /
							(  opt.mat[or,"UV.next"] - opt.mat[or,"UV"] - r[!exact] * (opt.mat[or,"L.next"]-opt.mat[or,"L"]))

  if (also.row.ind) {
    return(list(ind=ind,L=L))
  } else {
	  return(L)
  }
}            

