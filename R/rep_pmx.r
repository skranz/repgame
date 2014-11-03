# Search on a mesh
# We have an nx dimensional vector of activities x
# There exists functions c(x) and G(x) that generate cheating payoffs
# and best reply payoffs for vector of activity

# The nx dimensional vector grid.size stores the search grid size for every activity
# Basically, we will round activity x[j] to a multiple of grid.size[j]

# We have global ranges x.range for every activity

# A point is optimal for a given grid size if there is no other point
# that achieves a higher level of G (respectively a lower level of c[,i]) 
# with equal or lower liquidity

# Checks whether a set of new activities can improve on the optimal known activities

# opt.mat has the following structure: L, Y, within, x.... 

OUTSIDE  = 0
BOUNDARY = 1
WITHIN   = 2
WITHIN.OR.GLOBAL.BOUNDARY  = 3


pmx.get.opt.points = function(m,x.mat, opt.mat,k = "e", cnt, 
                              new.point.within.val = FALSE, tol.Y = REP.CNT$TOL) {
  restore.point("pmx.get.opt.points")
  # restore.point("pmx.get.opt.points")


  x.free = cnt$x.free;
  x.linked = cnt$x.linked;
  link.fun = cnt$link.fun;

  nx = m$nx;
  nx.free = NROW(x.free)
  
  # If some activities are linked to the values of other activities
  if (nx > nx.free) {
    if (NCOL(x.mat) == nx) {
      xf.mat = x.mat[,x.free,drop=FALSE]
    } else {
      xf.mat = x.mat
    }
    xl.mat = link.fun(xf.mat,k=k)
    if (NCOL(xl.mat) == m$nx) {
      x.mat = xl.mat
    } else {
      x.mat = matrix(NA,NROW(xf.mat),NCOL(xf.mat))
      x.mat[,x.free] = xf.mat
      x.mat[,x.linked] = xl.mat
    }
  }
  
  
  c.mat = m$cx.fun(x.mat)
  g.mat = m$gx.fun(x.mat)

  L.x = rowSums(c.mat-g.mat)
  
  ind  = which(L.x < -(10^-6))
  if (NROW(ind)>0) {
    mat = cbind(x.mat[ind,,drop=FALSE],L.x[ind],g.mat[ind,,drop=FALSE],c.mat[ind,,drop=FALSE])
    colnames(mat)=c(m$names.x,"L",paste("g",1:m$n,sep=""), paste("c",1:m$n,sep=""))
    print(mat)
    warning("negative L values found (see above), you probably misspecified cx.fun!")#
  }
    
  if (k=="e") {
    type = "Ue"
    sgn.Y = 1
    Y.x = rowSums(g.mat)
  } else {
    type = "vi"
    sgn.Y = -1
    Y.x = c.mat[,k]
  }
  
  if (cnt$use.L.lim) {
    rows = L.x >= cnt$L.min & L.x <= cnt$L.max
    L.x = L.x[rows]
    x.mat = x.mat[rows,]
    Y.x = Y.x[rows,]
  }
  
  if (!is.null(opt.mat)) {
    mat = rbind(cbind(L.x,Y.x,new.point.within.val,x.mat),
                 opt.mat)
    colnames(mat) = colnames(opt.mat)
    new.point = c(rep(TRUE,NROW(x.mat)),rep(FALSE,NROW(opt.mat)))
  } else {
    mat = cbind(L.x,Y.x,new.point.within.val,x.mat)
    new.point = rep(TRUE,NROW(x.mat))
  }
  
  ord = order(mat[,1],-sgn.Y*mat[,2],new.point,decreasing=FALSE)
  mat = mat[ord,,drop=FALSE]
  new.point = new.point[ord]
  
  if (sgn.Y == 1) {
    Y.max = c(-Inf,cummax(mat[-NROW(mat),2]))
    del = mat[,2] <= Y.max + tol.Y 
  } else {
    Y.min = c(Inf,cummin(mat[-NROW(mat),2]))
    del = mat[,2] >= Y.min - tol.Y
  }
  #cbind(mat,Y.max,del,mat[,2] - Y.max)[1:10,]
  
  return(list(opt.mat = mat[!del,,drop=FALSE], new.point = new.point[!del]))
  #return(mat[!del,,drop=FALSE])  
}

pmx.is.within.boundary.outside = function(x.mat, x.range, x.global.range) {
  
  # restore.point("pmx.is.within.boundary.outside")

  within.mat = boundary.mat = global.mat = matrix(NA,NROW(x.mat),NCOL(x.mat))
  for (j in 1:NCOL(x.mat)) {
    within.mat[,j]   = x.mat[,j] > x.range[j,1] & x.mat[,j] < x.range[j,2]
    boundary.mat[,j] = x.mat[,j] == x.range[j,1] | x.mat[,j] == x.range[j,2]
    global.mat[,j]  = x.mat[,j] <= x.global.range[j,1] | x.mat[,j] >= x.global.range[j,2]
  }
  within = rowSums(within.mat) == NCOL(x.mat)
  boundary = rowSums(within.mat + boundary.mat) == NCOL(x.mat) & (!within)
  within.or.global.boundary = rowSums(within.mat + global.mat) == NCOL(x.mat)
  
  ret = rep(OUTSIDE,NROW(x.mat))
  ret[boundary] = BOUNDARY
  ret[within.or.global.boundary] = WITHIN.OR.GLOBAL.BOUNDARY  
  ret[within]= WITHIN
  ret
}



pmx.get.within.rows = function(x.mat, x.range, x.global.range, old.within) {
  
  # restore.point("pmx.get.within.rows")

  all.within = old.within
  rows = which(all.within == FALSE)
  if (length(rows)>0) {
    mat = x.mat[rows,,drop=FALSE]
    within = rep(TRUE,NROW(mat))
    for (j in 1:NCOL(mat)) {
      within = within &  
               ((mat[,j] > x.range[j,1] & mat[,j] < x.range[j,2]) |
                (mat[,j] <= x.global.range[j,1] | mat[,j] >= x.global.range[j,2]))
    }
    all.within[rows] = within
  }
  return(all.within)
}


lseq = function(from,to,by) {
  lapply(1:NROW(from), function(i) seq(from[i],to[i],by[i]))
}


pmx.solve.game = function(m,cnt=NULL,cnt.e=cnt,cnt.i=cnt) {	           
	
	#restore.point("pm.solve.game")
	
	# remove old solution
	m$opt.mat = m$Ue.opt = m$vi.opt = NULL;
	
		
	ret = pmx.get.opt.mat(m,cnt.e=cnt.e,cnt.i=cnt.i)
	m$opt.mat = ret$opt.mat
	m$Ue.opt = ret$Ue.opt
	m$vi.opt = ret$vi.opt

	return(m)
}


refine.solution = function(m,cnt=NULL,cnt.e=cnt,cnt.i=cnt) {	           
	
	#restore.point("pm.refine.solution")
	
	if (m$sol.type != "pmx") {
	  warning("refine.solution so far only implemented for sol.type = pmx")
	  return(m)
  }
	
	ret = pmx.get.opt.mat(m,cnt.e=cnt.e,cnt.i=cnt.i)
	m$opt.mat = ret$opt.mat
	m$Ue.opt = ret$Ue.opt
	m$vi.opt = ret$vi.opt
	return(m)
}



pmx.get.opt.mat = function(m, symmetric = m$symmetric,
                        keep.only.opt.rows = !TRUE, do.plot=FALSE,cnt=NULL,cnt.e=cnt, cnt.i=cnt,...) {
	                        
	
	#restore.point("pmx.get.opt.mat")
	                        
	# Initialize matrices	                        
	if (symmetric) {
		i.max = 1
	} else {
		i.max = m$n
	}

	do.nothing = !is.null(cnt.e$do.nothing)
	if (do.nothing) do.nothing = cnt.e$do.nothing;	
	if (!do.nothing) {
	  Ue.opt = pmx.get.opt.k(m, k = "e",cnt=cnt.e)
  } else {
  	Ue.opt = m$Ue.opt  
  }
	min.L = max(0,min(Ue.opt[,"L"]))

	do.nothing = !is.null(cnt.i$do.nothing)
	if (do.nothing) do.nothing = cnt.i$do.nothing;	
	if (!do.nothing) {
  	vi.opt = list()
  	for (i in 1:i.max) {
  		#vi.opt[[i]] =  pmx.get.opt.k(m, k = i,...)
  	  vi.opt[[i]] =  pmx.get.opt.k(m, k = i,cnt=cnt.i)
  	  min.L = max(min.L, min(vi.opt[[i]][,"L"]))
  	}
	} else {
  	vi.opt = m$vi.opt
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
	# If the minimum L differs between states
	all.L = all.L[all.L>=min.L]
	
	all.L = sort(unique(all.L),decreasing=FALSE)
	Ue.opt.ind = findInterval(all.L,Ue.opt[,"L"])
	
	vi.opt.ind = matrix(NA,NROW(all.L),i.max)
	for (i in 1:i.max) {
		vi.opt.ind[,i] = findInterval(all.L,vi.opt[[i]][,"L"])
	}
	
	collab = c("delta","L","Ue","V",paste("v",1:m$n,sep=""),
	           paste("e.",m$names.x,sep=""),
	           paste(rep(1:i.max,each=m$nx),".", rep(m$names.x,times=i.max),sep=""),
	           "r","UV","opt")
	vi.start = which(collab=="v1"); vi.cols = paste("v",1:m$n,sep="")
	x.start =  vi.start+m$n; x.cols = x.start:(x.start+m$nx*i.max-1);
	mat = matrix(NA,NROW(all.L),NROW(collab))
	colnames(mat)=collab
	
	mat[,"L"] = all.L
	e.x.names = paste("e.",m$names.x,sep="")
	mat[,c("Ue",e.x.names)] = Ue.opt[Ue.opt.ind,c("Ue",m$names.x)]
	for (i in 1:i.max) {
  	i.x.names = paste(i,".",m$names.x,sep="")
		mat[,c(vi.cols[i],i.x.names)] = vi.opt[[i]][vi.opt.ind[,i],c("vi",m$names.x)]
	}

  
	if (i.max > 1) {
		mat[,"V"] = rowSums(mat[,vi.cols[1:i.max],drop=FALSE])
	} else {
		mat[,vi.cols[-1]] = mat[,"v1"]
		mat[,"V"] = m$n*mat[,"v1"]
	}
		
	mat[,"UV"] = (mat[,"Ue"]-mat[,"V"])
	mat[,"r"] = (mat[,"Ue"]-mat[,"V"]) / mat[,"L"]
	mat[!is.finite(mat[,"r"]),"r"] = Inf
	mat[,"delta"] = 1 / (1+mat[,"r"])

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
	
	list(opt.mat = mat, Ue.opt = Ue.opt, vi.opt = vi.opt)
}



# Initialize default values for refinement control
pmx.normalize.cnt = function(cnt,m,k) {
  
  # restore.point("pmx.normalize.cnt")
  
  nx = m$nx
  x.range = m$x.range
  n = m$n
  
  def.cnt = list(method="random", # "grid" is the alternative
                x.range = m$x.range, x.width = NULL,
                # For the random method
                size.draws = 50000, num.draws = 10,
                step.size = NULL,abs.dist = NULL,prob.draw.global = 0.2,
                local.rel.width = 0.05, local.abs.width = NULL,

                # For the adaptive grid refinement method
                step.size.start = NULL, step.size.end = NULL,
                num.step.start =    16, num.step.end = 1024,
                num.refinements = 1, step.change.factor = NULL,
                grid.size = 50000, start.grid.size = NULL, neighbor.grid.size = NULL,
                max.grid.size = 1000000,
                use.random.start.points = FALSE,num.random.start.points = 50000,
                dist.in.steps = NULL, min.dist.in.steps = 2,
                
                # General parameters
                do.nothing = FALSE, old.opt.within = 4, refine.opt.mat = TRUE,
                use.L.lim = FALSE, L.min = -Inf, L.max = Inf,
                # always.use.x.mat
                always.use.x.mat = NULL              
                )
   
  if (k=="e") {
    def.cnt$x.free = m$x.free.e;
    def.cnt$x.linked = m$x.linked.e;
    def.cnt$link.fun = m$link.fun.e;
  } else {
    k = as.numeric(k)
    def.cnt$x.free = m$x.free.i[[k]];
    def.cnt$x.linked = m$x.linked.i[[k]];
    def.cnt$link.fun = m$link.fun.i; 
  }

  
  # Copy defined parts of cnt into def.cnt
  def.cnt[names(cnt)] <- cnt               
  cnt <- def.cnt
  
  
  
  cnt$x.width = cnt$x.range[,2]-cnt$x.range[,1]

  set.to.default = function(parnames, def.value) {
    is.null = sapply(parnames, function(par) {
      if (is.null(cnt[[par]])){
        cnt[[par]] <<- def.value
      }})
  }
 
  set.to.default(c("start.grid.size","neighbor.grid.size"),cnt$grid.size)
  
  set.to.default(c("local.abs.width"),cnt$local.rel.width*cnt$x.width)
  set.to.default(c("abs.dist"),cnt$local.abs.width)

  set.to.default(c("step.size.start"),cnt$x.width / cnt$num.step.start)
  set.to.default(c("step.size.end"),cnt$x.width / cnt$num.step.end)

    
  cnt$num.step.start = round(cnt$x.width / cnt$step.size.start)
  cnt$num.step.end = round(cnt$x.width / cnt$step.size.end)
  
  # start* z^k = end <=> k = log(end/start,base=z)
  # z = (end / start)^(1/k)
  set.to.default(c("step.change.factor"), 
     (cnt$num.step.end / cnt$num.step.start)^(1/cnt$num.refinements))
  cnt$step.change.factor[cnt$step.change.factor <= 1] = 2
 
  assign.list(cnt)

  nx.free = NROW(cnt$x.free)
  cnt$step.size.start = rep(cnt$step.size.start, length.out = nx)
  cnt$step.size.end = rep(cnt$step.size.end, length.out = nx)
  cnt$step.change.factor = rep(cnt$step.change.factor, length.out = nx)
  cnt$step.size = rep(cnt$step.size, length.out = nx)
  
  if (is.null(cnt$dist.in.steps) & cnt$method=="grid") {
    # nx^(2*dist.in.steps+1) = grid.size
    # (2*dist.in.steps+1) = log_nx (grid.size)
    # 
    dist.in.steps = floor((log(cnt$neighbor.grid.size,base=nx.free)-1) / 2)
   
    if (dist.in.steps < cnt$min.dist.in.steps) {
      dist.in.steps = cnt$min.dist.in.steps
      str = (paste("dist.in.steps ",dist.in.steps," < cnt$min.dist.in.steps = ", cnt$min.dist.in.steps, " given ",
               " neighbor.grid.size ", cnt$neighbor.grid.size, ". Grid.size with min.dist.in.steps = ",
               nx.free^(2*dist.in.steps+1)))
      warning(str)
      print(str)
      flush.console()
    }
    cnt$dist.in.steps = dist.in.steps
  }    
  cnt$dist.in.steps = rep(cnt$dist.in.steps, length.out = nx)
  cnt
}



pmx.get.opt.k = function(m, k = "e", cnt= NULL) { 
  restore.point("pmx.get.opt.k")
  # restore.point("pmx.get.opt.k")
  
  cnt = pmx.normalize.cnt(cnt,m,k)  
  assign.list(cnt)  # Copies all elements of cnt into the local environment
 
  opt.mat = pmx.get.opt.mat.from.m(m=m,k=k,cnt=cnt) 
  
  if (cnt$method == "random") {
    return(pmx.get.opt.k.random(m=m,k=k,cnt=cnt))
  }
    
  if (k=="e") {
    Y.lab = "Ue"
  } else {
    Y.lab = paste("vi")
  }
  nx = m$nx;
        
  
  step.size = step.size.start
  x.grid = lseq(x.range[,1],x.range[,2],by = step.size)
  if (!use.random.start.points) { 
    xf.mat = make.grid.matrix(x=x.grid[x.free])
  } else {
    xf.mat = pmx.draw.x.mat(num.random.start.points, x.grid = x.grid[x.free], replace=TRUE)
  }
  if (!is.null(cnt$always.use.xmat)) {
    xf.mat = rbind(cnt$always.use.xmat,xf.mat)
  }
  
  ret = pmx.get.opt.points(m=m,x.mat = xf.mat,cnt=cnt, opt.mat=opt.mat,k = k)
  opt.mat = ret$opt.mat

  colnames(opt.mat) = c("L",Y.lab,"within",m$names.x)

  if (!use.random.start.points) {
    step.size.above = step.size > step.size.end
    if (sum(step.size.above)==0)
      return(opt.mat)
    step.size = step.size[step.size.above] / step.change.factor[step.size.above]
  }
 
  plot(opt.mat[,1],opt.mat[,2], main=paste("State",k," Step size=", round(step.size[1],5)),
         xlab="L",ylab="G(ae) / c.i(ai)")
  
  counter = 0 
  while (TRUE) {
    
    step.size.above = step.size > step.size.end - 10^(-9)
    if (sum(step.size.above[x.free]) == 0)
      break;
  
    opt.mat[,"within"] = FALSE
    
    
    while (sum(!opt.mat[,"within"]) > 0) {
      counter = counter+1
      if (cnt$use.L.lim) {
        rows = which((!opt.mat[,"within"]) &
                     opt.mat[,"L"] >= cnt$L.min & opt.mat[,"L"] <= cnt$L.max)
      } else {
        rows = which(!(opt.mat[,"within"]))
      }
      if (length(rows)==0)
        break();
        
      #ind = sample(rows,1)
      ind = rows[1]
      neigh = pmx.get.neighborhood(opt.mat[ind,-c(1:3)][x.free],
             step.size[x.free],x.range[x.free,,drop=FALSE],dist.in.steps = dist.in.steps[x.free])
      ret = pmx.get.opt.points(m=m,x.mat = neigh$x.mat,cnt=cnt, opt.mat=opt.mat,k=k)
      opt.mat = ret$opt.mat
      
      #place = pmx.is.within.boundary.outside(opt.mat[,-c(1:3),drop=FALSE][,x.free,drop=FALSE],
      #                                   neigh$x.range, x.range[x.free,,drop=FALSE])
      #opt.mat[,"within"] = opt.mat[,"within"] | place == WITHIN |
      #                     place == WITHIN.OR.GLOBAL.BOUNDARY 
      
      opt.mat[,"within"] = pmx.get.within.rows(opt.mat[,-c(1:3),drop=FALSE][,x.free,drop=FALSE],
                                         neigh$x.range, x.range[x.free,,drop=FALSE], opt.mat[,"within"])

      points(opt.mat[ret$new.point,1],opt.mat[ret$new.point,2],col="blue")
      if (counter %% 3 == 0)
        try(legend("right",legend=paste(counter," NROW(opt.mat)=",NROW(opt.mat)),
                   bg="lightgrey"))

       
    }
    plot(opt.mat[,1],opt.mat[,2], main=paste("State",k," Step size=", round(step.size[1],5)),
         xlab="L",ylab="G(ae) / c.i(ai)")
    
    step.size[step.size.above] = step.size[step.size.above] / step.change.factor[step.size.above]
 
  }
  opt.mat
} 

pmx.get.neighborhood = function(x,step.size,x.range, dist.in.steps = rep(3,NROW(x))) {
  restore.point("pmx.get.neighborhood")
  # restore.point("pmx.get.neighborhood")
  
  x.min = pmax(x - step.size * dist.in.steps, x.range[,1])
  x.max = pmin(x + step.size * dist.in.steps, x.range[,2])
  
  x.grid = lseq(x.min,x.max,by = step.size)
  x.mat = make.grid.matrix(x=x.grid)
  return(list(x.mat=x.mat, x.range = cbind(x.min,x.max)))
}

pmx.get.opt.mat.from.m = function(m,k,cnt) {
  restore.point("pmx.get.opt.mat.from.m")
  # restore.point("pmx.get.opt.mat.from.m")
  if (!cnt$refine.opt.mat) 
    return(NULL)
  
  k.names.x = paste(k,".",m$names.x,sep="")
  if (k=="e") {
    return(m$Ue.opt)
  } else {
    k = as.numeric(k)
    if (is.null(m$vi.opt))
      return(NULL)
    return(m$vi.opt[[k]])
  }
}

pmx.get.opt.k.random = function(m, k = "e", cnt) { 
  restore.point("pmx.get.opt.k.pure.random")
  # restore.point("pmx.get.opt.k.pure.random")
 
  opt.mat = pmx.get.opt.mat.from.m(m,k,cnt)
  
  if (cnt$do.nothing) 
    return(opt.mat)

  nx = m$nx; 
  assign.list(cnt)
  if (k=="e") {
    Y.lab = "Ue"
  } else {
    k = as.numeric(k)
    Y.lab = paste("vi")
  }
  global.x.range = x.range

  last.skipped = FALSE
  for (d in 1:num.draws) {
    draw.global = (runif(1) < prob.draw.global | is.null(opt.mat) | last.skipped)
    last.skipped = FALSE
    if (draw.global) {
      xf.mat = pmx.draw.x.mat(size.draw, x.range = global.x.range[x.free,], replace=TRUE)
      new.point.within.val = 0
    } else {
      min.within = min(opt.mat[,"within"])
      if (!is.finite(min.within)) {
        last.skipped = TRUE
        next();
      }
      if (cnt$use.L.lim) {
        rows = which(opt.mat[,"within"]==min.within &
                     opt.mat[,"L"] >= cnt$L.min & opt.mat[,"L"] <= cnt$L.max)
      } else {
        rows = which(opt.mat[,"within"]==min.within)
      }
      if (NROW(rows)==0) {
        last.skipped = TRUE
        next();
      }
      ind = sample(rows,1)
      x = opt.mat[ind,-c(1,2,3)]
      opt.mat[ind,"within"] = Inf
      new.point.within.val = min.within+1
      x.min = pmax(x - abs.dist, global.x.range[,1])
      x.max = pmin(x + abs.dist, global.x.range[,2])
      xf.mat = pmx.draw.x.mat(size.draw,  x.range = cbind(x.min,x.max)[x.free,], replace=TRUE)
    }
    if (!is.null(cnt$step.size))
      xf.mat = round.to.step(xf.mat,cnt$step.size[1])
      
    if (d == 1 & !is.null(cnt$always.use.xmat))
      xf.mat = rbind(cnt$always.use.xmat,xf.mat)

    ret = pmx.get.opt.points(m=m,x.mat = xf.mat, opt.mat=opt.mat,k = k, cnt=cnt,new.point.within.val = new.point.within.val)
    opt.mat = ret$opt.mat
    if (d==1) {
      colnames(opt.mat) = c("L",Y.lab,"within",m$names.x)
    }
    if (d %% 10 == 1 | d <6) {
      ylim = range(opt.mat[,2])
      ylim[2] = ylim[2]+ 0.25*(ylim[2]-ylim[1])
      plot(opt.mat[,1],opt.mat[,2], main=paste("State",k),
         xlab="L",ylab="G(ae) / c.i(ai)", ylim=ylim)
    } else {
      points(opt.mat[,1],opt.mat[,2])
    }
    try(legend("topright",legend=paste("Draws ",d," n.opt",NROW(opt.mat)),
                   bg="lightgrey"))
  }
  opt.mat
} 


#' Inits a game of perfect monitoring with large or infinite action space
#'
#' See the tutorial "Interactively Solving Repeated Games" for examples
pmx.init.game = function(n=2,gx.fun=NULL,cx.fun=NULL, cix.fun=NULL, Gx.fun = NULL, Lx.fun = NULL,
    nx = NROW(x.range),x.range, names.x = NULL, ix=NULL, Xi = NULL, x.free = NULL, link.fun = NULL,
    x.free.e = x.free, x.free.i = x.free, 
    link.fun.e = link.fun, link.fun.i = link.fun,
    name=NULL, symmetric=FALSE) {
      
  
	#restore.point("pmx.init.game")
	
	repgames.startup()
	
	m = list()	  
  class(m)=c("pmx","repgame","list")
  m$name = name
  m$sol.type = "pmx"
  m$n = n
  m$nx = nx
	m$symmetric = symmetric
	m$i.max = m$n
	if (symmetric)
	  m$i.max = 1
	
  if (is.null(names.x))
	  names.x = paste("x",1:nx,sep="")
  m$names.x = names.x
  
  if (!is.matrix(x.range)) {
    x.range = matrix(x.range,nx,2,byrow=TRUE)
  } else if (NROW(x.range)==1) {
    x.range = matrix(as.vector(x.range),nx,2,byrow=TRUE)
  }
  colnames(x.range)=c("min.x","max.x")
  rownames(x.range)=m$names.x
  m$x.range = x.range
  
  ########################################################################
  # Store info about linked variables
  ########################################################################
  
  if (is.null(x.free.e)) {
    m$x.free.e = 1:nx
  } else {
    m$x.free.e = x.free.e
  }
  m$x.linked.e = setdiff(1:nx,m$x.free.e)
  m$nx.free.e = NROW(m$x.free.e)
  m$nx.linked.e = NROW(m$x.linked.e)
  m$link.fun.e = link.fun.e
  
  m$link.fun.i = link.fun.i
  if (is.null(x.free.i)) {
    m$x.free.i = list()
    for (i in 1:m$i.max) {
      m$x.free.i[[i]] = m$x.free.e
    }
  } else {
    if (is.list(x.free.i)) {
      m$x.free.i  = x.free.i
    } else {
      m$x.free.i = list()
      for (i in 1:m$i.max) {
        m$x.free.i[[i]] = x.free.i
      }      
    }
  }
  m$x.linked.i = list()
  m$nx.free.i = numeric(m$n)
  m$nx.linked.i = numeric(m$n)
  for (i in 1:m$i.max) {
    m$x.linked.i[[i]] = setdiff(1:nx,m$x.free.i[[i]])
    m$nx.free.i[i] = NROW(m$x.free.i[[i]])
    m$nx.linked.i[i] = NROW(m$x.linked.i[[i]])
  }
  
  ################################################################
  # Mapping of activities to players
  ################################################################
  
#   if (is.null(ix) | is.null(ix)) {
#     if (nx==n) {
#       ix = 1:m$n
#     } else {
#       stop("No. of activities x unequal to number of players. You must specify Xi or ix to match activities to players")
#     }
#   }
#   
#   if (is.null(ix)) {
#     ix = rep(NA,nx)
#     for (i in 1:m$n) 
#       ix[Xi[[i]]] = i
#   }
#   m$ix = ix
#   names(m$ix) = m$x.names
# 
#   if (is.null(Xi)) {
#     Xi = list()
#     for (i in 1:m$n) 
#       Xi[[i]] = which(ix == i)
#   }
#   m$Xi = Xi
	  
  # Payoff functions and cheating payoff functions
    
  m$gx.fun = gx.fun
  m$cx.fun = cx.fun
    
	if (is.null(Gx.fun))
  	Gx.fun = function(x.mat) rowSums(gx.fun(x.mat))
  m$Gx.fun = Gx.fun
  	
  if (is.null(Lx.fun))
  	Lx.fun = function(x.mat) rowSums(cx.fun(x.mat)-gx.fun(x.mat))
	m$Lx.fun = Lx.fun
	                             
  return(m)
}

# Generate randome ranges
# Width, and start point (or end point) of the range are
# randomly determined
generate.local.range = function(num,min=0,max=1) {
  width = runif(num,min,max*0.7)
  inner = runif(num,min,max)
  left = sample(c(TRUE,FALSE),num,replace=TRUE)
  right=!left
  range = matrix(NA,num,2)
  range[left,1] = pmax(0,inner-width)[left]
  range[left,2] = pmax(inner,width)[left]
  range[right,1] = pmin(1-width,inner)[right]
  range[right,2] = pmin(inner+width,1)[right]
  range
}
generate.local.range(3)
  

#method can take "unif", "two.range" or "sample"

pmx.draw.x.mat = function(ndraw, x.range=NULL,x.grid = NULL, step.size=NULL  , replace=TRUE, method=NULL) {
  
  # restore.point("pmx.draw.x.mat")
  
  if (is.null(method)) {
    method = "sample"
    if (is.null(x.grid) & is.null(step.size))
      method = "two.range"
  }
  
  # Two range is the default procedure. It allows positive or negative correlation
  # patterns between activities  
  if (method=="two.range") {
    nx = NROW(x.range) 
    x.mat  = matrix(NA,ndraw,nx)
    # Draw two ranges with start and endpoints between 0 and 1
    loc.range = generate.local.range(2)
    # Define how many activities shall lie in range 1 (the others lie in range 2)
    num.range1 = sample(1:nx,1)
    # Determine which of the activities lie in range 1
    which.range = sample(c(rep(1,num.range1),rep(2,nx-num.range1)),nx)
    
    # Map the ranges between 0 and 1 to the actual ranges of each activity
    x.range[,1] = x.range[,1]+loc.range[which.range,1]*(x.range[,2]-x.range[,1])
    x.range[,2] = x.range[,1]+loc.range[which.range,2]*(x.range[,2]-x.range[,1])
    
    # Draw activities from uniform distribution between the ranges determined above
    for (j in 1:nx) {
      x.mat[,j]= runif(ndraw,x.range[j,1],x.range[j,2])
    }
    return(x.mat)    
  } else if (method=="sample") {
    if (is.null(x.grid)) {
      if (is.null(step.size)) {
        step.size = x.range[,2]-x.range[,1] / 10000
      }
      x.grid = lseq(x.range[,1],x.range[,2],by = step.size)
    }
       
    nx = length(x.grid)
    x.mat  = matrix(NA,ndraw,nx)
    for (j in 1:nx) {
      x.mat[,j] = sample(x.grid[[j]],ndraw,replace = replace)
    }
    return(x.mat)
  } else {
    nx = NROW(x.range)
    x.mat  = matrix(NA,ndraw,nx)
    for (j in 1:nx) {
      x.mat[,j]= runif(ndraw,x.range[j,1],x.range[j,2])
    }
    return(x.mat)
  }
}

round.to.step = function(x,step) {
  round(x / step) * step
}  

#' Numerically test whether a supplied stage game best reply function seems to be correct
#'
#' See the tutorial for details
check.cx.fun = function(m, num = 50, i.of.x = 1:m$n,k="e", x.sample = "random",
                        method="grid", num.grid.steps=1000, use.i = 1:m$n, i.vec = NULL) {
  
  # restore.point("check.cx.fun")
  
  if (NROW(i.of.x) != m$nx)
    stop(paste("Model has ",m$nx, " activities, you need to specify a vector i.of.x of this size!"))
  if (k=="e") {
    x.free = m$x.free.e
  } else {
    x.free = m$x.free.i
  }
  is.free = 1:m$nx %in% x.free
  
  x.mat = matrix(NA,num,m$nx);
  if (x.sample=="random") {
    for (p in 1:num) {
      x.mat[p,] = pmx.draw.x.mat(1, m$x.range)
    }
  } else if (x.sample=="opt.mat") {
    cols = paste(k,".",m$names.x,sep="")
    x.mat = m$opt.mat[,cols,drop=FALSE]
  } else {
    x.mat = x.sample
  }
  num = NROW(x.mat)
  
  if (is.null(i.vec)) {
    i.sample = intersect(unique(i.of.x[x.free]),use.i)
    i.vec = sample(i.sample,num,replace=TRUE)
  }
  i.vec = rep(i.vec,length.out=NROW(x.mat))
  
  
  mat = matrix(NA,num,4+2*m$nx)
  colnames(mat)=c("cx","cn","diff","i",paste("org",m$names.x),paste("br",m$names.x))
  cn = cx = rep(NA,num)
 
  for (p in 1:num) {
    i = i.vec[p] 
    x = x.mat[p,,drop=FALSE]
    f = function(x) {
      if (!is.matrix(x)) {
        x = t(x)
        colnames(x) = m$names.x
      }
      m$gx.fun(x)[,i]
    }
    fp = which(i.of.x == i & is.free)
    ret = sk.optim(as.vector(x),f,lower = m$x.range[,1], upper = m$x.range[,2],
                   free.par = fp,num.grid.steps=num.grid.steps, method=method,
                   f.can.take.matrix = TRUE)
    mat[p,"cn"] = ret$value
    mat[p,"cx"] = m$cx.fun(x)[,i]
    mat[p,"i"] = i
    mat[p,5:(4+m$nx)] = as.vector(x)
    mat[p,(5+m$nx):(4+2*m$nx)] = ret$par
  }    
  mat[,"diff"] = mat[,"cx"]-mat[,"cn"]
  if (x.sample=="random") {
    mat=mat[order(mat[,"cn"]),]
  }
  ylim = range(mat[,c("cn","cx")])
  plot(1:num,mat[,"cx"],col="black",ylim=ylim,ylab="",xlab="",main="Cheating Payoffs: Analytical vs Numerical")
  points(1:num,mat[,"cn"],col="blue",pch=3)
  return(mat)
}


make.numerical.cx.fun = function(m,i.of.x = 1:m$n,method="grid",
                                 num.grid.steps=1000,gx.needs.colnames=FALSE,tol=0.0001) {
  restore.point("cx.fun")
  cx.fun = function(x.mat) {
    restore.point("cx.fun")
     #restore.point("cx.fun");restore.point("make.numerical.cx.fun")
    c.mat = matrix(NA,NROW(x.mat),m$n)
    for (p in 1:NROW(x.mat)) {
      for (i in 1:m$n) {
        x = x.mat[p,,drop=FALSE]
        f = function(x) {
          if (!is.matrix(x)) {
            x = t(x)
            if (gx.needs.colnames)
              colnames(x) = colnames(x.mat)
          }
          m$gx.fun(x)[,i]
        }
        fp = which(i.of.x == i)
        ret = sk.optim(as.vector(x),f,lower = m$x.range[,1], upper = m$x.range[,2],
                       free.par = fp,num.grid.steps=num.grid.steps, method=method,
                       f.can.take.matrix = TRUE,tol=tol)
        c.mat[p,i] = ret$value
      }
    }
    return(c.mat)
  } 
  return(cx.fun)
}
   
   
#     
#     
# pmx.check.cx.fun = function(m, num.points = 10, step.size = NULL) {
#   if (is.null(step.size)) {
#     step.size = (m$x.range[,2]-m$x.range[,1]) / 1000
#   }
#   x.mat = pmx.draw.x.mat(num.points, m$x.range, step.size)
#   i.test = sample(1:m$n,NROW(x.mat))
#   cx.from.fun = m$cx.fun(x.mat)
#   
#   cx.from.test = matrix(NA,NROW(x.mat),m$n)
#   
#   f.fun = function(xi,i,x) {
#     x[i] = xi
#     x = as.matrix(x,1,NROW(x))
#     m$gx.fun(x)[,i]
#   }
#   
#   for (r in 1:NROW() {
#     cx.from.text = optimize(f.fun,interval = m$x.range[
#   


