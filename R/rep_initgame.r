
###################################################################################################################
#  Functions related to creation of a model m
###################################################################################################################
			
make.cheating.payoffs = function(m,add.labels = TRUE) {
  
	#restore.point("make.cheating.payoffs")
		
	# Method for 2 player games
	if (m$n==2 & !is.null(m$g2)) {	
		#c1 = apply(m$g1,2,max)
		c1 = colMaxs(m$g1) 
		if (m$symmetric) {
			c2 = c1
		} else {
  		c2 = rowMaxs(m$g2)
			#c2 = apply(m$g2,1,max)
		}
		c = cbind(rep(c1,times=NROW(m$g1)),rep(c2,each=NCOL(m$g1)))
		return(c)
	}
		
  m$action.gi = GridInd(dim=m$a.dim)
  c.mat = matrix(NA,m$nA,m$n)  
  for (i in 1:m$n) {
    repl = make.replies.gi(m$action.gi,k.i = i)
    g.replies = matrix(m$g[repl$reply.mat,i],nrow=NROW(repl$reply.mat))
    g.max = rowMaxs(g.replies)
    c.mat[,i] = g.max[repl$cp.ind]
  }  
  return(c.mat)

	
	
	# Method for n player games, considerably slower than the two player method. Even though it is coded in C
	if(!is.loaded("repgames"))
		library.dynam("repgames", package = "repgames", lib.loc = NULL,
              verbose = getOption("verbose"),
              file.ext = .Platform$dynlib.ext)


		
	m$c = rep(-99999,m$nA*m$n)	

	ret = .C("make_cheating_payoff",as.integer(m$n), as.integer(m$a.dim), as.integer(m$shift.ai),
	  as.double(m$g),as.double(m$c))

	c = matrix(ret[[5]],m$nA,m$n)

	plot(c.mat-c)
		
	return(c)  
}

#' Init a new repeated game, see the tutorial for a description and examples
init.game = function(n=2,g=NULL,ai.dim = NULL, g.fun=NULL,action.val=NULL, g1=NULL, g2 = NULL, c.fun = NULL,
										 phi.mat=NULL,name="", symmetric=FALSE, lab.ai = NULL,
                     lab.a = NULL, lab.y =NULL, add.labels = TRUE, audit.prob = 0,
                     add.prob.to.1=FALSE, sol.type = NULL) {
	  
		#restore.point("init.game")
	
		repgames.startup()
		  
		m = list()	  
	  class(m)=c("repgame","list")
	  m$name = name
	  if (is.null(sol.type)) {
  	  if (is.null(phi.mat)) {
    	  sol.type="pm"
  	  } else {
    	  sol.type = "imp"
  	  }
	  }
	  m$sol.type = sol.type
	  m$n = n
	  
  	                     		
	  action.val.list = NULL
		if (!is.null(g) & !is.null(ai.dim)) {
			gtype="g.mat"
		} else if (!is.null(g.fun) & !is.null(action.val)) {                 	                             
			gtype="g.fun"
			if (!is.list(action.val)) {
				action.val.list = replicate(m$n,action.val, simplify=FALSE)
			} else {
				action.val.list = action.val
			}
		} else if (!is.null(g1)) {
			gtype="g1"
			if (is.null(g2)) {
				g2 = t(g1)
				symmetric=TRUE
			}
		} else {
			stop("init.game: You have to specify either g and ai.dim or g.fun and action.val or g1.")
		}

			
	                             
	  m$symmetric = symmetric
		m$audit.prob = audit.prob
	  
	  # Initialize ai.dim
	  if (gtype=="g.mat")  {
		  if (NROW(ai.dim)==1) {
				m$ai.dim = rep(ai.dim,m$n)
			} else {
				m$ai.dim = ai.dim
			}
		} else if (gtype=="g.fun") {
			if (is.list(action.val)) {
				m$a.dim = sapply(action.val,NROW)
			} else {
				m$a.dim = rep(NROW(action.val),m$n)
			}
		} else if (gtype=="g1") {
			m$a.dim = c(NROW(g1),NCOL(g1))
		}

		
		m$shift.ai = rev(cumprod(rev(m$a.dim))) / m$a.dim[1:n]
		m$modulus.ai = c(prod(m$a.dim),m$shift.ai[-n])
		# Number of action profiles
		m$nA = prod(m$a.dim)
		
		# Initialize labels of the different actions
		if (!is.null(lab.ai)) {
			if (is.list(lab.ai)) {
				m$lab.ai = lab.ai
				if (length(lab.ai)==1)
					for (i in 2:m$n) m$lab.ai[[i]] = m$lab.ai[[1]];
				
			} else {
				m$lab.ai = replicate(m$n,lab.ai,simplify=FALSE)
			}
			if (NROW(m$lab.ai[[1]]) != m$a.dim[1]) {
				stop("Error: the parameter lab.ai has the wrong dimensions.")
			}
			if (gtype=="g1") {
				rownames(g1)=rownames(g2)=m$lab.ai[[1]]
				colnames(g1)=colnames(g2)=m$lab.ai[[2]]			
			}
		} else {
			m$lab.ai = list()
		  if (gtype=="g.mat")  {
			  # If gtype=="g.mat" then lab.ai must be supplied
			} else if (gtype=="g.fun") {
				if (!is.null(names(action.val.list[[1]]))) {
					for (i in 1:n) m$lab.ai[[i]] = names(action.val.list[[i]]);
				} else {
					m$lab.ai = action.val.list
				}
			} else if (gtype=="g1") {
				m$lab.ai[[1]] = rownames(g1)
				m$lab.ai[[2]] = colnames(g1)
			}
		}

		
		
		# Initialize labels of the different action profiles
		if (add.labels) { 
			if (is.null(lab.a)) {
				if (n==2) {
					# Default format. In a prisoners' dilemma game looks e.g. "C|D"
	  			m$lab.a = as.vector(t(outer(m$lab.ai[[1]],m$lab.ai[[2]],function(x,y){paste(x,y,sep="|")})))
  			} else {
	  			lab.a.work = make.grid.matrix(x=m$lab.ai)
	  			m$lab.a = paste(lab.a.work[,1],lab.a.work[,2],sep="|")
	  			for (i in 3:m$n)
	  				m$lab.a = paste(m$lab.a,lab.a.work[,i],sep="|")
  			}
			} else {
				m$lab.a = lab.a
			}
			# m$ind.a can be used to find the index of an action profile if its name is known
			m$ind.a = 1:m$nA
			names(m$ind.a) = m$lab.a
		}
		
		# A vector that specifies whether a certain action profile is symmetric or not		
		m$sym.a = rep(FALSE,m$nA)
		if (n ==2) {
			if(m$a.dim[1] == m$a.dim[2]) {
				m$sym.a[((1:m$a.dim[1])-1)*m$a.dim[1]+(1:m$a.dim[1])]=TRUE
			}
		}
		####################################################################################################
		# Generate payoffs								
		####################################################################################################
		
		if (is.null(action.val.list)) 
  	  action.val.list = lapply(m$a.dim, function(x) 1:x)
		
  	m$action.val.mat = make.grid.matrix(x=action.val.list)
  	
		
		# Initialize ai.dim
	  if (gtype=="g.mat")  {
		  m$g = g
		} else if (gtype=="g.fun") {
			val.mat = make.grid.matrix(x=action.val.list)
			m$g = g.fun(val.mat)
		} else if (gtype=="g1") {
			m$g1 = g1
			m$g2 = g2
			# Table that contains all payoffs nA rows and n columns
			m$g = cbind(as.vector(t(g1)),as.vector(t(g2)))
		}

		if (add.labels) {
			rownames(m$g)=m$lab.a
		}
		m$G = rowSums(m$g)
		
		# Generate cheating payoffs
		if (!is.null(c.fun) & !is.null(g.fun)) {
			m$c = c.fun(val.mat)
		} else {
			m$c = make.cheating.payoffs(m)
		}
		if (add.labels) {
			rownames(m$c)=rownames(m$g)
		}

		m$C = rowSums(m$c)

		# Vector of liquidity requirements. Will be partially calculated if the model is solved (for imperfect monitoring)
		m$L = rep(NA,m$nA);
		#if (only.per.mon) {
		#	m$L = m$C-m$G
		#}
		
		
		m$nash = which(m$C == m$G)
		m$v.nash = rep(NA,m$n)
		if (NROW(m$nash)>0)
			for (i in 1:m$n) m$v.nash[i] = min(m$g[m$nash,i])
		
		m$v.min.max = apply(m$c,2,min)
		m$V.min.max = sum(m$v.min.max)
		m$L[m$nash] = 0		
		
		
		if (add.labels) {
			names(m$sym.a) = names(m$L) = rownames(m$g)
		}		
			
		if (!is.null(phi.mat)) {
	  	m$multi.phi.mat = FALSE
			if (is.vector(phi.mat)) {
				phi.mat = matrix(phi.mat,1,length(phi.mat))
			}
			sum.prob = apply(phi.mat,2,sum)		
			if (add.prob.to.1 & sum(sum.prob < 1) > 0 ) {
				str = rownames(phi.mat)
				phi.mat = rbind(phi.mat,1-sum.prob)
				rownames(phi.mat) = c(str,"yL")
			}
			m$ny = NROW(phi.mat)
			m$phi.mat = phi.mat
			m$lab.y = lab.y
			if (is.null(lab.y)) {
				if (!is.null(rownames(phi.mat))) {
					m$lab.y = rownames(phi.mat)
				} else {
					m$lab.y = paste("y",1:m$ny,sep="")
				}
			}
			rownames(m$phi.mat) = m$lab.y
			colnames(m$phi.mat)= m$lab.a
			
		}

	  m$lp=m$lp.info = m$lp.i = m$lambda = NULL
			
		# Init default lp where all liquidity is given to player n
		if (!is.null(phi.mat)) {
			ret = make.default.lp.for.m(m,i=m$n)
			m$lp = ret$lp
			m$lp.info = ret$lp.info
			m$lp.i = ret$lp.info$i
			m$lambda = ret$lp.info$lambda
		}
		return(m)
}

	

get.dummy.phi.mat = function(ny,nA) {
	mat = matrix(0,ny,nA)
	mat[1,]=1
	return(mat)
}

#' Returns a model that uses grim.trigger punishments instead of optimal penal codes
set.to.grim.trigger = function(m) {
  mg = m
  if (mg$sol.type == "pmx" | mg$sol.type == "imp") {
    warning("Can only set models of sol.type pm grim.trigger")
    return(m)
  }
  mg$sol.type = "grim"
  mg$Ue.opt = mg$vi.opt = NULL
  mg$opt.mat = NULL
  mg$phi.mat = NULL
  mg
}

#' Returns a model that uses Abreu's symmetric stick and carrot punishments as in Abreu (1986) instead of optimal (possibly asymmetric) penal codes with transfers
set.to.Abreu.stick.carrot = function(m) {
  ma = m
  if (ma$sol.type == "pmx" | ma$sol.type == "imp") {
    stop("Can only set models of sol.type pm to Abreu.stick and carrot")
  }
  if (!m$symmetric) {
    stop("The game must be symmetric (m$symmetric==TRUE)")
  }

  ma$sol.type = "Abreu"
  ma$Ue.opt = ma$vi.opt = NULL
  ma$opt.mat = NULL
  ma$phi.mat = NULL
  ma
}

#' Set signal distribution of a model with imperfect public monitoring
set.phi.mat = function(m,phi.mat, lab.y = rownames(phi.mat)) {
  rownames(phi.mat) = lab.y
  colnames(phi.mat) = m$lab.a
  m$ny = NROW(phi.mat)
  m$phi.mat = phi.mat
  m$lab.y = rownames(phi.mat)
  
  
  ret = make.default.lp.for.m(m,i=m$n)
	m$lp = ret$lp
	m$lp.info = ret$lp.info
	m$lp.i = ret$lp.info$i
	m$lambda = ret$lp.info$lambda

	m$sol.type = "imp"
  m
}

#' See the Tutorial
set.audit.prob = function(m, audit.prob) {
  m$audit.prob = audit.prob
  m
}

make.new.lp.for.m = function(m) {
  ret = make.default.lp.for.m(m,i=m$n)
	m$lp = ret$lp
	m$lp.info = ret$lp.info
	m$lp.i = ret$lp.info$i
	m$lambda = ret$lp.info$lambda
  m
}

add.phi.mat.li = function(m,phi.mat.li, lab.y = NULL) {
	m$multi.phi.mat = FALSE
	phi.mat = phi.mat.li[[1]]
	if (is.null(lab.y)) {
		lab.y = rownames(phi.mat.li[[1]])
	}
	m$lab.y = lab.y
	m$ny = NROW(phi.mat)
	m$phi.mat.li = phi.mat.li
	m$sol.type = "imp"
	m
}

# Let A be a list with the names of all action profiles
make.phi.mat = function(A,trembles=NULL) {
	# If A is payoff matrix with rownames and colnames
	if (is.matrix(A)) {
		A = as.vector(t(outer(rownames(A),colnames(A),function(x,y){paste(x,y,sep="")})))
	}
	# If A is just a number
	if(length(A)==1 & is.numeric(A)) {
		A = 1:A
		rownames(A) = 1:length(A)
	}
	nA = length(A)
	phi.mat = diag(nA)
	colnames(phi.mat)=A
	rownames(phi.mat)=paste("y",A,sep="")
	if (is.null(trembles)) {
		return(phi.mat)
	}
	if (!is.numeric(trembles[,1])) {
		tr = matrix(0,NROW(trembles),3)
		tr[,1] = match(trembles[,1],A)
		tr[,2] = match(trembles[,2],A)
		tr[,3] = as.numeric(trembles[,3])
	} else {
		tr = trembles;
	}
	colnames(tr)=c("source","dest","prob")
		
	for (t in 1:NROW(tr)) {
		phi.mat[tr[t,1],tr[t,1]] = phi.mat[tr[t,1],tr[t,1]] - tr[t,3]
		phi.mat[tr[t,2],tr[t,1]] = phi.mat[tr[t,2],tr[t,1]] + tr[t,3]
	}
	if (max(phi.mat)>1 | min(phi.mat)<0) {
		warning("Entries in phi.mat that are no probabilities!")
	}
	zero.row = (apply(phi.mat,1,sum) == 0)
	phi.mat = phi.mat[!zero.row,]
	return(phi.mat)
}

