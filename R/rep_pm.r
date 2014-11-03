
###########################################################################################################
# Finds for games with perfect monitoring the list of optimal equilibrium state profiles and corresponding
# optimal action profiles as function of the total liquidity L
###########################################################################################################


##################################################
# An unelegant function based on pmin
##################################################
rowMins = function(x) {
  # Construct a call pmin(x[,1],x[,2],...x[,NCOL(x)])
	code = paste("x[,",1:(NCOL(x)),"]",sep="",collapse=",")
	code = paste("pmin(",code,",na.rm=TRUE)")
	return(eval(parse(text=code)))
}



# A quick function to return minimum of each row
# The code is unelegant but much quicker than calling
# apply(x,1,max)
rowMaxs = function(x) {
  # Need to construct a call pmax(x[,1],x[,2],...x[,NCOL(x)])
	code = paste("x[,",1:(NCOL(x)),"]",sep="",collapse=",")
	code = paste("pmax(",code,",na.rm=TRUE)")
	return(eval(parse(text=code)))
}


get.grim.trigger.mat = function(m) {
	
	#restore.point("get.grim.trigger.mat")

	if (NROW(m$nash)<1) {
		warning("get.grim.trigger.mat: the stage game has no Nash equilibrium")
		return(NULL)
	}
	
	# Liquidity requirement
	L.i = m$c-m$g
	# Neccessary discount rate and factor for player i
	r.i = L.i
	for (i in 1:m$n) {
		r.i[,i] = (m$g[,i] - m$v.nash[i]) / L.i[,i]
	}
	r.i[m$nash,] = Inf
	delta.i = 1 / (1+r.i)
	
	#r = apply(r.i,1,min)
	r  = rowMins(r.i)
	#delta = apply(delta.i,1,max)
	delta = rowMaxs(delta.i)
	#L = apply(L.i,1,max)
	L = rowMaxs(L.i)
	
	
	A.ord = order(m$G,m$sym.a,decreasing=TRUE)
	mat = cbind(delta[A.ord],m$G[A.ord],sum(m$v.nash),r[A.ord],(1:m$nA)[A.ord],
	            delta.i[A.ord,],m$g[A.ord,],L.i[A.ord,],L,1)
	colnames(mat) = c("delta","Ue","V","r","ae",
	                   paste("delta",1:m$n,sep=""),paste("g",1:m$n,sep=""),
	                   paste("L",1:m$n,sep=""),"L","opt")
	rownames(mat) = m$lab.a[A.ord]
	del = mat[,"r"] <= 0
	mat = mat[!del,,drop=FALSE]

	cum.min.delta = c(Inf,cummin(mat[,"delta"])[-NROW(mat)])
	del = mat[,"delta"] >= cum.min.delta
	mat = mat[!del,,drop=FALSE]
	
	
	
	return(mat[NROW(mat):1,,drop=FALSE])
}

# Get optimal 
pm.get.Ue.max = function(m, ignore.a = NULL) {
	
	#restore.point("pm.get.Ue.max")
	
	if (is.null(ignore.a)) {
  	max.a.check = m$nA
  	if (m$symmetric) {
  		A.ord = order(m$G,-(m$C-m$G),m$sym.a,decreasing=TRUE)
  	} else {
  		A.ord = order(m$G,-(m$C-m$G),decreasing=TRUE)
  	}
	} else {
  	max.a.check = sum(!ignore.a)
  	if (m$symmetric) {
  		A.ord = order(-ignore.a,m$G,-(m$C-m$G),m$sym.a,decreasing=TRUE)
  	} else {
  		A.ord = order(-ignore.a,m$G,-(m$C-m$G),decreasing=TRUE)
  	}
  	A.ord = A.ord[1:max.a.check]
	}
  		
	cumL = c(Inf,cummin(m$C[A.ord]-m$G[A.ord])[-NROW(A.ord)])
	A.inc = c(m$C[A.ord]-m$G[A.ord] < cumL)
	A.ord = A.ord[A.inc]	
	LU.PM = cbind(m$C[A.ord]-m$G[A.ord],m$G[A.ord],A.ord)	
	collab = c("L","Ue","a")
	colnames(LU.PM) = collab
	rownames(LU.PM) = rownames(m$g)[A.ord]

	#plot(LU.PM)
	return(LU.PM[NROW(LU.PM):1,,drop=FALSE])
}

# Get optimal punishment profiles
pm.get.vi.min = function(m,i=1, ignore.a = NULL) {
  if (is.null(ignore.a)) {
    max.a.check = m$nA
  	A.ord = order(-m$c[,i],-(m$C-m$G),decreasing=TRUE)
	} else {
  	max.a.check = sum(!ignore.a)
  	A.ord = order(-ignore.a,-m$c[,i],-(m$C-m$G),decreasing=TRUE)
  	A.ord = A.ord[1:max.a.check]
	}

	cumL = c(Inf,cummin(m$C[A.ord]-m$G[A.ord])[-NROW(A.ord)])
	A.inc = c(m$C[A.ord]-m$G[A.ord] < cumL)
	A.ord = A.ord[A.inc]	
	Lvi = cbind(m$C[A.ord]-m$G[A.ord],m$c[A.ord,i],A.ord)	
	collab = c("L","vi","a")
	colnames(Lvi) = collab
	rownames(Lvi) = rownames(m$g)[A.ord]
	return(Lvi[NROW(Lvi):1,,drop=FALSE])
}	

pm.solve.game = function(m, keep.only.opt.rows = FALSE, ignore.ae = NULL, ignore.a1 = NULL,ignore.ai = NULL) {	           
	
	#restore.point("pm.solve.game")
	
	
	if (m$sol.type=="grim") {
  	m$opt.mat = get.grim.trigger.mat(m)
	} else if (m$sol.type=="Abreu") {
  	stopifnot(m$symmetric)
  	if (!is.null(ignore.ae)) {
  	  ignore.ae =  (!m$sym.a) & ignore.ae
	  } else {
  	  ignore.ae =  (!m$sym.a)
	  }
  	if (!is.null(ignore.a1)) {
  	  ignore.a1 =  (!m$sym.a) & ignore.a1
	  } else {
  	  ignore.a1 =  (!m$sym.a)
	  }
 	  m$opt.mat = pm.get.opt.mat(m, keep.only.opt.rows = keep.only.opt.rows,
	                  ignore.ae = ignore.ae, ignore.a1 = ignore.a1,ignore.ai = ignore.ai,
	                  Abreu.stick.carrot=TRUE)

  } else {
	  m$opt.mat = pm.get.opt.mat(m, keep.only.opt.rows = keep.only.opt.rows,
	                  ignore.ae = ignore.ae, ignore.a1 = ignore.a1,ignore.ai = ignore.ai)
  }
	return(m)
}

# Adds end points of horizontal lines, so that we can simply the resulting matrix
pm.make.plot.mat = function(opt.mat) {
	
	#restore.point("pm.make.plot.mat")

	
	plot.mat = opt.mat[rep(1:NROW(opt.mat),each=2),,drop=FALSE]
	rows = (1:NROW(opt.mat))*2
	#plot.mat[rows,"opt"] = 0  # Does not work for grim.trigger.mat
	plot.mat[rows[-NROW(rows)],c("L","delta","r")]=plot.mat[rows[-NROW(rows)]+1,c("L","delta","r")]
	plot.mat[NROW(plot.mat),c("delta","r")]=c(1,0)
	return(plot.mat)	
}

pm.get.opt.mat = function(m,Ue.opt = NULL, vi.opt =NULL, symmetric = m$symmetric,
                        keep.only.opt.rows = TRUE, do.plot=FALSE, 
                        ignore.ae = NULL, ignore.a1 = NULL,ignore.ai = NULL, Abreu.stick.carrot = FALSE) {
	                        
	
	#restore.point("pm.get.opt.mat")
	                        
	# Initialize matrices	                        
	if (symmetric) {
		i.max = 1
	} else {
		i.max = m$n
	}

	if (is.null(Ue.opt)) {
		Ue.opt = pm.get.Ue.max(m, ignore.a = ignore.ae)
	}

	if (is.null(vi.opt)) {
		for (i in 1:i.max) {
  		if (!is.null(ignore.ai)) {
    		vi.opt[[i]] =  pm.get.vi.min(m,i=i, ignore.a = ignore.ai[[i]])
      } else {
        vi.opt[[i]] =  pm.get.vi.min(m,i=i, ignore.a = ignore.a1)
      }

		}
	}
	
	
	# Initialize vectors that merge information from optimal matrices of all states.
	n.L = NROW(Ue.opt)
	for (i in 1:i.max) n.L = n.L + NROW(vi.opt[[i]]);
	all.L = rep(NA,n.L);
	rows = 1:NROW(Ue.opt)
	all.L[rows] = Ue.opt[,"L"];
	for (i in 1:i.max) {
		rows = max(rows)+1:NROW(vi.opt[[i]])
		all.L[rows] = vi.opt[[i]][,"L"];
	}
	
	all.L = sort(unique(all.L),decreasing=FALSE)
	Ue.opt.ind = findInterval(all.L,Ue.opt[,"L"])
	
	vi.opt.ind = matrix(NA,NROW(all.L),i.max)
	for (i in 1:i.max) {
		vi.opt.ind[,i] = findInterval(all.L,vi.opt[[i]][,"L"])
	}
	
	collab = c("delta","L","Ue","V",paste("v",1:m$n,sep=""),"ae",paste("a",1:m$n,sep=""),
	           "r","UV","opt")
	vi.start = which(collab=="v1"); vi.cols = vi.start:(vi.start+m$n-1);
	ai.start = which(collab=="a1"); ai.cols = ai.start:(ai.start+m$n-1);
	mat = matrix(NA,NROW(all.L),NROW(collab))
	colnames(mat)=collab
	
	mat[,"L"] = all.L
	mat[,c("Ue","ae")] = Ue.opt[Ue.opt.ind,c("Ue","a")]
	for (i in 1:i.max) {
		mat[,c(vi.cols[i],ai.cols[i])] = vi.opt[[i]][vi.opt.ind[,i],c("vi","a")]
	}
	if (i.max > 1) {
		mat[,"V"] = rowSums(mat[,vi.cols[1:i.max],drop=FALSE])
	} else {
		mat[,vi.cols[-1]] = mat[,"v1"]
		mat[,"V"] = m$n*mat[,"v1"]
	}
		
	mat[,"UV"] = (mat[,"Ue"]-mat[,"V"])
	if (!Abreu.stick.carrot) {
	  mat[,"r"] = (mat[,"Ue"]-mat[,"V"]) / mat[,"L"]
	  mat[!is.finite(mat[,"r"]),"r"] = Inf
	  mat[,"delta"] = 1 / (1+mat[,"r"])
  } else if (Abreu.stick.carrot) {
  	# Works only for strongly symmetric equilibria
  	# No monetary transfers are conducted
  	
  	# Let V.sc = (1-delta) * mat[,"V"] + delta * mat[,"Ue"]
  	# We need
  	# Ue >= (1-delta) * U.Br(Ue) + delta * V.sc
  	# Ue >= (1-delta) * (Ue + L) + delta * ((1-delta)V + delta*Ue)
  	# (delta - delta^2) Ue >= (1-delta)L + delta*(1-delta)*V
  	# (1-delta)*delta *Ue >= (1-delta)L + delta*(1-delta)*V
  	# delta*Ue >= L + delta*V
  	# delta >= L / (Ue-V)
  	
  	delta = mat[,"L"] / (mat[,"Ue"]-mat[,"V"])
		delta[mat[,"L"] == 0] = 0
		delta[delta>1] = 1
		mat[,"delta"] = delta
		mat[,"r"] = (1-delta) / delta
	}


	# The lowest discount factor of all action structures with strictly higher L
	min.delta.right = rev(cummin(c(1,rev(mat[,"delta"]))))[-1]
	# Only keep those elements that have a strictly lower critical discount factor than elements with higher L
	mat[,"opt"]= mat[,"delta"]<min.delta.right-REP.CNT$TOL.PM
	#cbind(min.delta.right,mat)
	if (keep.only.opt.rows) {
		mat = mat[mat[,"opt"]==1,,drop=FALSE]
	}
	if (!is.null(m$lab.a)) {
		lab = paste("(",m$lab.a[mat[,"ae"]],")",sep="")
		for (i in 1:i.max) {
			lab = paste(lab,",(",m$lab.a[mat[,ai.cols[i]]],")",sep="")
		}
		rownames(mat) = lab
	}
	
	mat
}



init.quantity.competition.PM = function(n,q.seq,P.fun,c.fun, P.hom = !is.list(P.fun), P.sym = P.hom,name="Quantity Competition") {
	q.li = NULL; P.li = NULL; c.li = NULL;
	q.sym = TRUE; c.sym = TRUE; P.hom = TRUE;
	if (is.list(q.seq)) {
		q.li = q.seq
		a.dim = sapply(q.li,length)
		a.li = lapply(a.dim, function(d) {1:d})
		q.sym = FALSE;
	} else {
		a.dim = rep(NROW(q.seq),n)
	}
	if (is.list(P.fun)) {P.li = P.fun;P.hom = FALSE;}
	if (is.list(c.fun)) {c.li = c.fun;c.sym = FALSE;}
	
	sym = q.sym & P.sym & c.sym
	
	if (q.sym) {
		am = qm = make.grid.matrix(n,1:NROW(q.seq))
		for (i in 1:n) {
			qm[,i] = q.seq[am[,i]]
		}
	} else {
		am = qm = make.grid.matrix(n,a.li)
		for (i in 1:n) {
			qm[,i] = q.li[[i]][am[,i]]
		}
	}
	nA = NROW(am)
	
	g = matrix(NA,nA,n)
	if (P.hom) {
		Q = rowSums(qm)
		P = P.fun(Q)
	}
	for (i in 1:n) {
		if (!P.hom) {P.fun = P.li[[i]]}
		if (!c.sym) {c.fun = c.li[[i]]}			
		if (!P.hom) {
			P = apply(qm,1,P.fun)
		}
		g[,i] = P*qm[,i]-c.fun(qm[,i])
	}
	lab = round(qm[,1],2)
	for (i in 2:n) {lab = paste(lab,round(qm[,i],2),sep="|")}
	rownames(g)=lab
	m = list(name=name,perfect.monitoring = TRUE,symmetric=sym, n = n,a.dim = a.dim, qm = qm, g=g)
	m = init.model.PM(m)
}	

# Homogeneous Cournot
# nq = 501;
# n = 2; q.seq = seq(0,100 / (n-1),length=nq);
# P.fun=function(Q) {pmax(100-Q,0)};
# c.fun = function(qi){5*qi};
# P.sym = TRUE; name = "Quant. Comp"
# 	
# m = init.quantity.competition.PM(n=n,q.seq=q.seq,P.fun=P.fun,c.fun=c.fun, P.hom=TRUE,
#                                  P.sym = TRUE,name="Quantity Competition")

pm.get.opt.mat.row = function(m,delta) {
	max(which(m$opt.mat[,"delta"]<=delta))
}

# This function delivers an optimal payment plan for a given
# the construction follows the construction of optimal payments
# in the proof of Proposition ### in Kranz & Ohlendorf (2010)
pm.get.optimal.payment.plan = function(m, delta, ap=NULL) {
	
	#restore.point("pm.get.optimal.payment.plan")

	if (is.null(ap)) {
		row = pm.get.opt.mat.row(m,delta)
		collab = colnames(m$opt.mat)
		vi.start = which(collab=="v1"); vi.cols = vi.start:(vi.start+m$n-1);
	  ai.start = which(collab=="a1"); ai.cols = ai.start:(ai.start+m$n-1);
		ap = m$opt.mat[row,c(ai.start-1,ai.cols)]
	}
	names(ap)=c("ae",paste("a",1:m$n,sep=""))
	ae = ap[1];ai = ap[-1];
	
	L.k = m$C[ap]-m$G[ap]
	lambda.k = matrix(NA,m$n+1,m$n)
	rownames(lambda.k) = names(ap)
	for (k in 1:(m$n+1)) {
		lambda.k[k,] = (m$c[ap[k],]-m$g[ap[k],]) / L.k[k]
	}
	lambda.k[L.k==0,] = rep(1/m$n,m$n)
	
	# Punishment payoffs
	v = m$c[cbind(ai,1:m$n)]
	
	# L.star and lambda.star
	L.star = (delta/(1-delta))*(m$G[ae]-sum(v))
	lambda.star = delta* (lambda.k["ae",] + ((m$g[ae,] - v) / L.star))

	# Payments if complying or deviating
	p.comp= matrix(NA,m$n+1,m$n)
	rownames(p.comp)=names(ap)
	p.dev = p.comp


	for (k in 1:(m$n+1)) {
		p.comp[k,] = 0 + (lambda.star - lambda.k[k,]) * L.star
		p.dev[k,] = p.comp[k,] + m$c[ap[k],] - m$g[ap[k],]
	}

# 	ue             = m$g[ae,]-p.comp["ae",]
# 	up.front.equal = (ue - (sum(ue)/m$n)) / (1-delta)
	
	names(ap) = m$lab.a[ap]
	ret = list(
		delta = delta,
		ap = ap,
		p.dev = p.dev,
		p.comp = p.comp
# 		,ge = m$g[ae,],
# 		ce = m$c[ae,],
#    	up.front.equal = up.front.equal,
# 		ue = m$g[ae,]-p.comp["ae",]

	)
	return(ret)	
}


	
                                        
# Make defection signal
make.cheating.phi.mat.li = function(m,ind.err.prob=0.1,ind.detect.prob=1/2,
                                      ano.err.prob=0.1,ano.detect.prob=1/2,ind.signal=TRUE, ano.signal=TRUE ) {
	if (!ind.signal) {
		ind.err.prob = 0;ind.detect.prob = 0;
	}
	if (!ano.signal) {
		ano.err.prob = 0;ano.detect.prob = 0;
	}
	
	if (ind.detect.prob + ano.detect.prob > 1) {
		warning("make.cheating.phi.mat.li: ind.detect.prob + ano.detect.prob shall not exceed 1. Automatically scaled.")
		ind.detect.prob = ind.detect.prob / (ind.detect.prob + ano.detect.prob)
		ano.detect.prob = ano.detect.prob / (ind.detect.prob + ano.detect.prob)
	}
		                                      
  y.ok = rep(1-ind.detect.prob-ano.detect.prob,m$nA)
  y.ano.cheat = rep(ano.detect.prob,m$nA)
  y.ind.cheat = matrix(0,m$n,m$nA)
  
  phi.mat = rbind(y.ok,y.ano.cheat,y.ind.cheat)
  rownames(phi.mat) = c("yok","yano",paste("y",1:n,sep=""))
  colnames(phi.mat) = rownames(m$g)
  
  org.phi.mat = phi.mat
  phi.mat.li = vector("list",m$nA)
  a = 1
  
  # Note that the signals after multilateral deviations do not add up to 1
  # But that is irrelevant, since we do not have incentive constraints for 
  # multilateral deviations
  for (a in 1:m$nA) {
	  # y.ok[a] = 1-ind.err.prob-ano.err.prob
	  phi.mat[1,a] = 1-ind.err.prob-ano.err.prob
	  #y.ano.cheat[a] = ano.err.prob
	  phi.mat[2,a] = ano.err.prob
	  phi.mat[3:NROW(phi.mat),a]=ind.err.prob / m$n
	  
	  for (i in 1:n) {
		  a.ind = a
		  a.rep = setdiff(get.replies.ind(m,a.ind=a,i=i),a)
			#y.ind.cheat[i,a.rep] = ind.detect.prob / m$n
			phi.mat[i+2,a.rep] = ind.detect.prob
		}
		phi.mat.li[[a]] = phi.mat
		phi.mat = org.phi.mat
	}
	names(phi.mat.li) = rownames(m$g)
	m$phi.mat.li = phi.mat.li
	m$ny = NROW(phi.mat)
	m$lab.y = rownames(phi.mat)
	return(m)
}
