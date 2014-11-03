# Copies as list into an environment
assign.list = function(li, envir=sys.frame(-1),clone.env=TRUE) {
  if (is.null(li)) {return()}
	
	if (!clone.env) {
  	for (i in 1:length(li)) {
  		assign(names(li)[i],li[[i]],envir=envir)
  	}
	} else {
  	for (i in 1:length(li)) {
    	obj = li[[i]]
    	if (is.environment(obj)) {
      	obj = clone.environment(obj)	
		  }
		  assign(names(li)[i],obj,envir=envir)
  	}
  }	
}


matrix.get.rows = function(mat,col.val,col.names=names(col.val))  {
 
  restore.point("matrix.get.rows")
  
  if (!is.matrix(mat))
    mat = as.matrix(mat)
  if (!is.list(col.val))
    col.val = as.list(col.val) 
  mwhich(mat,cols=col.names,vals=col.val,'eq','AND')
}


plot.multi.lines = function(mat=NULL,xvar,yvar,ynames=yvar,col=NULL,ylim=NULL,xlab=xvar,ylab="",
                               legend.pos=NULL,legend.title=NULL,add=FALSE,lwd=1,...) {
  
  if (is.null(ylim)) {
		ylim = range(mat[,yvar])
  	if (!is.null(legend.pos)) {
      #if (legend.pos=="topleft" | legend.pos == "topright") {
          #ylim[2] = ylim[2]+0.2*diff(ylim)
      #}
    }               		
	}
  ny = NROW(yvar)
  if (is.null(col)) {
    if (ny<=5) {
      col=c("blue","red","green","black","orange")[1:ny]
    } else {
      col = rainbow(ny)
    }
	}

	if (!add) {
  	plot(mat[,xvar],mat[,yvar[1]],type="n",ylim=ylim,xlab=xlab,ylab=ylab,...)
	}
	# Draw lines
 	for (i in 1:NROW(yvar))
		lines(mat[,xvar],mat[,yvar[i]], col=col[i],lwd=lwd)

	# Draw lines once more dotted, so that we can better see lines that are on top of each other better
	if (NROW(yvar)>1) {
 		# Draw lines
 		for (i in (NROW(yvar)-1):1) 
			lines(mat[,xvar],mat[,yvar[i]], col=col[i],lty=2,lwd=lwd)
	}

	# Draw legend, if desired	
	if (!is.null(legend.pos)) {
		legend(legend.pos, legend=ynames, fill=col,title=legend.title)
	}
}


# Get a matrix with all rows that are unique in cols[1] and add a count
get.unique.table = function(dat,cols, add.count = TRUE) {
  un = unique(dat[,cols[1]])
  rows = match(un,dat[,cols[1]])
  mat = cbind(dat[rows,cols],as.numeric(table(dat[,cols[1]])))
  colnames(mat) = c(colnames(dat[1:2,cols]),"count")
  mat
}

# Generate a matrix of polynomials of order o
cross.poly.mat = function(mat,cols=colnames(mat),ord=2) {
  nc = NROW(cols)
  li = replicate(nc,list(0:ord))
  grid = expand.grid(li)
  fun.pol = function(pol) {
    ret = 1
    for (k in 1:NROW(pol))
      ret = ret * mat[,k]^pol[k]
    ret
  }
  poly = apply(grid,1,fun.pol)
  name.fun = function(pol,names=1:NROW(pol)) {
    use = pol>0
    paste(names[use],"^",pol[use],sep="",collapse="")      
  }
  colnames(poly)=apply(grid,1,name.fun,names=cols)
  poly = poly[,-1]
  poly
}

# Call for debugging
#checkUsageEnv(globalenv())

# Looks through all loaded functions and searches for
# global variables that are used within the functions
# this is a common source for errors
check.global.vars = function() {
  require(codetools)
  print("Usage of global variables in loaded functions:")
  funs = ls.funs(env=globalenv())
  for (fn in funs) {
    cmd = paste("findGlobals(",fn,",merge=FALSE)$variables",sep="")
    glob = try(eval(parse(text=cmd)))
    if (length(glob)>0) {
      print("")
      print(paste(paste(fn,": "),paste(glob,collapse=", ")))
    }
  }
}
#check.global.vars()


with.floor = function(mat,floor=0) {
  mat[mat<floor] = floor
  mat
}

with.floor.and.ceiling = function(mat,floor,ceiling) {
  mat[mat<floor] = floor
  mat[mat>ceiling] = ceiling
  mat
}


with.ceiling = function(mat,ceiling) {
  mat[mat>ceiling] = ceiling
  mat
}



sublist = function(li, cols, simplify = TRUE) { 
  sl = list()
  if (!simplify | length(cols) > 1) {
    for (i in 1:length(li)) {
      sl[[i]] = li[[i]][cols]
    }
  } else {
    for (i in 1:length(li)) {
      sl[[i]] = li[[i]][[cols]]
    }
  }    
  return(sl)
}
# a = list(x=1:10,y=c("Hi","you!"),z=0)
# b = a
# li = list(a,b)
# sublist(li,cols=c("x","z"))
# sublist(li,"y")  


rbind.list = function(li, cols=NULL, only.common.cols=FALSE) {
  if (length(li)==1) {
    return(li[[1]])
  }
  if (only.common.cols) {
    names = colnames(li[[1]])
    for (i in 2:length(li)) {
      names = intersect(names,colnames(li[[i]]))
    }
    if (length(names)<1) {
      warning("Matrices in list have no common columns")
      return(NULL)
    }
    len = sapply(li,NROW,simplify=TRUE)
    mat = matrix(NA,sum(len),length(names))
    colnames(mat) = names
    rowstart = 1
    for (i in 1:length(li)) {
      rowend = rowstart+len[i]-1
      mat[rowstart:rowend,] = li[[i]][,names]
      rowstart = rowend+1
    }
    return(mat)
  }
  
  if (is.null(cols)) {
    return(do.call("rbind",li))
  } else {
    return(do.call("rbind",li)[,cols])
  }
}

cbind.list = function(li) {
    return(do.call("cbind",li))
}
x = list(1:10,11:20)
rbind.list(x)



sk.findInterval = function(x, vec, match.smaller=TRUE, vec.ordered = "increasing",...) {

  #restore.local.objects("sk.findInterval")
  
  if (vec.ordered == "increasing") {
    ind = findInterval(x,vec,...)
    if (!match.smaller)
      ind = ind +1
    return(ind)
  } else {
    vec.ord = order(vec)
    #ind.ord = findInterval(x,vec[vec.ord],...)
    ind.ord = findInterval(x,vec[vec.ord])
   
    if (!match.smaller) {
      same = x %in% vec
      ind.ord[!same] = ind.ord[!same] +1
    }
   
    return(vec.ord[ind.ord])
  }
}

#sk.findInterval(1:10,vec = c(0,5,1,12,9), match.smaller = FALSE, vec.ordered = FALSE)


paste.matrix.cols = function(mat,cols=1:NCOL(mat),...) {
  if (NROW(cols)==2) {
    return(paste(mat[,cols[1]],mat[,cols[2]],...))
  } else if (NROW(cols)==3) {
    return(paste(mat[,cols[1]],mat[,cols[2]],mat[,cols[3]],...))
  } else {
    code = paste("mat[,",cols,"]",collapse=",")
    code = paste("paste(",code,",...)",sep="")
    return(eval(parse(text=code)))
  }
}

paste.matrix.rows = function(mat,rows=1:NROW(mat),...) {
  if (NROW(rows)==2) {
    return(paste(mat[rows[1],],mat[rows[2],],...))
  } else if (NROW(rows)==3) {
    return(paste(mat[rows[1],],mat[rows[2],],mat[rows[3],],...))
  } else {
    code = paste("mat[",rows,",]",collapse=",")
    code = paste("paste(",code,",...)",sep="")
    return(eval(parse(text=code)))
  }
}


# mat = matrix(1:100,10,10)
# paste.matrix.cols(mat)
# paste.matrix.rows(mat,sep=",")


add.rowvec = function(m,v) {
  return(t(t(m) +v))
} 

#' APPROXEQ Are a and b approximately equal (to within a specified tolerance)?
#' p = approxeq(a, b, thresh)
#' 'tol' defaults to 1e-3.
approxeq = function(a, b, tol=1e-3) {
  isTRUE(all.equal(a,b,tol=tol, check.attributes=FALSE))
}
all.eq = function(...) {
  isTRUE(all.equal(...,check.attributes=FALSE))
}



# Code seems more complicated than neccessary, but I tried to avoid
# parent.env(cenv)<- ...
# as they may become deprecated 
clone.environment = function(env, made.clones=ListEnv(org = list(), copy = list()), clone.parents=TRUE, clone.global = FALSE, exclude = NULL, clone.children = TRUE) {

	penv = parent.env(env)
  #Clone parents
  if (clone.parents & (!(identical(penv,emptyenv())
                    | (identical(penv,globalenv()) & ! clone.global)))) {
    cpenv = clone.environment(penv, made.clones = made.clones, clone.parents = TRUE)
  } else {
    cpenv <- penv
  } 
  
  if (length(made.clones$org) > 0) {
    for (i in 1:length(made.clones$org)) {
      if (identical(made.clones$org[[i]],env)) {
        return(made.clones$copy[[i]])
      }
    }
  }
  
  cenv = new.env(parent=cpenv)
  
  ind = length(made.clones$org)+1
  made.clones$org[[ind]] = env
  made.clones$copy[[ind]] = cenv
  
  #Clone children
  names  = setdiff(ls(env),exclude)
	#browser()
  if (clone.children) {
    for (na in names) {
      
			# If an error occurs here, the variable [[
			obj = env[[na]]
			
			
      if (is.environment(obj)) {
        obj = clone.environment(obj, made.clones = made.clones, clone.parents = TRUE, clone.global = clone.global)
        cenv[[na]] <- obj
      } else { 
        cenv[[na]] <- obj
      }
    }
  } else {
    for (na in names) {
      cenv[[na]] <- env[[na]]
    }
  }    
  class(cenv) <- class(env)
  return(cenv)  
}

clone.env = function(...) clone.environment(...)
# 
# env1 = ListEnv(a=10, b="Hi")
# env2 = new.ListEnv(parent=env1)
# env2$c = "c"
# env1$env = env2
# env2$a
# cenv = clone.environment(env2)

ls.funs <-function(env=sys.frame(-1)) {
  unlist(lapply(ls(env=env),function(x)if (is.function(get(x)))x))
}
ls.vars <- function(env=sys.frame(-1)) {
  unlist(lapply(ls(env=env),function(x)if(!is.function(get(x)))x))
}

currentenv = function() {
  sys.parent(1)
}



copy.env = function(dest=sys.frame(sys.parent(1)),source=sys.frame(sys.parent(1)),
        names = NULL, name.change = NULL, exclude=NULL) {
  #store.local.objects()  # DO NOT !!!!
  #restore.local.objects("copy.env")
  
  if (is.null(name.change)) {
    if (is.environment(source)) {
      if (is.null(names))
        names = ls(envir=source)
      names = setdiff(names,exclude)
      for (na in names) {
        assign(na,get(na,envir=source), envir=dest)
      }
    }
    if (is.list(source)) {
      if (is.null(names))
        names = names(source)
      names = setdiff(names,exclude)
      for (na in names) {
        assign(na,source[[na]], envir=dest)
      }
    }
  } else {
    if (is.environment(source)) {
      if (is.null(names))
        names = ls(envir=source)
      names = setdiff(names,exclude)
      for (na in names) {
        ind = match(na,name.change[,1])
        if (!is.na(ind)) {
          assign(name.change[ind,2],get(na,envir=source), envir=dest)
        } else {         
          assign(na,get(na,envir=source), envir=dest)
        }
      }
    }
    if (is.list(source)) {
      if (is.null(names))
        names = names(source)
      names = setdiff(names,exclude)
      for (na in names) {
        ind = match(na,name.change[,1])
        if (!is.na(ind)) {
          assign(name.change[ind,2],source[[na]], envir=dest)
        } else {         
          assign(na,source[[na]], envir=dest)
        }
      }
    }
  }      
}

copy.into.env = function(dest=sys.frame(sys.parent(1)), source=sys.frame(sys.parent(1)), names = NULL, name.change = NULL,exclude=NULL) {
 copy.env(dest,source,names,name.change,exclude)
}

copy.into.list = function(dest=NULL, source=sys.frame(sys.parent(1)), names = NULL,exclude=NULL,overwrite=TRUE) {
  if (is.null(dest)) {
    if (is.list(source)) {
      return(source)
    } else {
      dest=list()
    }
  }
  stopifnot(is.list(dest))

  if (!overwrite) {
    excluder = c(exclude,names(dest))
  }
  if (is.environment(source)) {
    if (is.null(names))
      names = ls(envir=source)

    names = setdiff(names,exclude)
    for (na in names) {
      dest[[na]]=get(na,envir=source)
    }
  }
  if (is.list(source)) {
    if (is.null(names))
      names = names(source)
    names = setdiff(names,exclude)
    for (na in names) {
      dest[[na]]=source[[na]]
    }
  }
  return(dest)
}


set.default = function(env,name,x, overwrite.null = TRUE, inherits = TRUE) {
  if (is.environment(env)) {
    if (!exists(name,envir=env, inherits = inherits)) {
      env[[name]] <- x
    } else {
      if (overwrite.null & is.null(env[[name]])) {
        env[[name]] <- x
      }
    }
    return(env)
  } else if (is.list(env)) {
    if (overwrite.null) {
      if (is.null(env[[name]])) {
        env[[name]] <- x
      }
    } else if (!(name %in% names(env))) {
      env[[name]] <- x
    }
    return(env)
  }
}
  

# Calculates the 2dimensional paretofrontier of the points val1 and val2
# The function returns the indices of the points that lie on the Pareto Frontier
# ordered by val1 and val2.
sk.pareto.frontier = function(val1,val2, tol=0, ord = NULL) {
  if (is.null(ord))
    ord = order(val1,val2,decreasing=TRUE)
  
  val2 = val2[ord]
	cummax2 = c(-Inf,cummax(val2[-NROW(val2)]))
	val2.inc = val2 > cummax2 + tol
	ord = ord[val2.inc]	
  return(ord)
}



# Finds position where the function f becomes zero
# First tries find.root and if this fails tries optimize
findzero = function(f, lower, upper, tol = .Machine$double.eps*10,result.tol = tol, try.uniroot=TRUE,...) {
  
  if (try.uniroot) {
    ret = tryCatch(uniroot(f,lower=lower,upper=upper,...), error = function(e) NULL)
    if (!is.null(ret)) {
      return(ret$root)
    }
  }

  f.sqr = function(...) f(...)^2
    
  ret = tryCatch(optimize(f.sqr,lower=lower,upper=upper,tol=tol,...), error = function(e) NULL)
  if (is.null(ret)) {
    warning("findzero: error in optimize")
    return(NA)
  } else if (abs(ret$objective)>result.tol) {
    warning("findzero: no solution found with optimize (min = ",ret$objective," > result.tol=",result.tol,")")
    return(NA)
  }
  ret$min
}

# A wrapper for optimization. Allows to specify which variables shall be free
# Has the same syntax for one and multidimensional optmization
# Uses optim, omptimize or a grid search
sk.optim = function(par,f, lower = NULL, upper = NULL, free.par = 1:NROW(par), method="default",
                    num.grid.steps = NULL,maximize=TRUE,f.can.take.matrix = FALSE,tol=.Machine$double.eps^0.25,...) {

  #restore.point(""sk.optim")
  
  
  n.par = NROW(par)
  if (!is.numeric(free.par)) 
    free.par = which(free.par)
  n.free = NROW(free.par)
  if (n.free == 0) {
    warning("No free parameters!")
    return(list(par=par,value=f(par,...)))
  }
  if (n.free != n.par) {
    sgn = 1
    if (maximize & n.free > 1 & method != "grid")
      sgn = -1
    g = function(x,...) {
      if (is.matrix(x)) {
        mat = matrix(par,NROW(x),NROW(par),byrow=TRUE)
        mat[,free.par] = x
        return(sgn*f(mat,...))
      } else {
        p = par
        p[free.par]=x
        return(sgn*f(p,...))
      }
    }
  } else {
    g = f
  }
  
  if (method=="default") {
    if (free.par == 1 & !is.null(lower)) {
      method = "optimize"
    } else if (!is.null(lower)) {
      method = "L-BFGS-B"
    } else {
      method = "Nelder-Mead"
    }
  }
      
      
  
  if (method == "grid") {
    if (is.null(num.grid.steps))
      stop("You have to specify num.grid.steps if you are using the grid method")
    num.grid.steps = rep(num.grid.steps,length.out=n.par)
    steps.li = lapply(free.par, function(i) seq(lower[i],upper[i],length=num.grid.steps[i]))
    
    par.mat = make.grid.matrix(x=steps.li)
    if (f.can.take.matrix) {
      val = g(par.mat,...)
    } else {
      val = sapply(1:NROW(par.mat), function(ind) g(par.mat[ind,],...))
    }
    if (maximize) {
      opt.ind = which.max(val)
    } else {
      opt.ind = which.min(val)
    }
    par.opt = par
    par.opt[free.par]= par.mat[opt.ind,]
    return(list(par=par.opt, value = val[opt.ind]))
  }
 
  
  # Use optimize
  if (n.free == 1) {
    ret = optimize(g,lower = lower[free.par], upper = upper[free.par],
                  maximum = maximize,tol=tol,...)
    par.opt = par 
    if (maximize) {
      par.opt[free.par] = ret$maximum
    } else { 
      par.opt[free.par] = ret$minimum
    }
    return(list(opt.ret = ret, par=par.opt, value = ret$objective))
  }
  
  if (method == "L-BFGS-B") {
    ret = optim(par[free.par],g,method=method,tol=tol,...)
  } else {
    if (is.null(lower)) {
      ret = optim(par[free.par],g,method=method,tol=tol)
    } else {
      stop("Only methods grid and L-BFGS-B so far implemented for constrained, multivariable optimization")
    }
  } 
  par.opt = par 
  par.opt[free.par] = ret$par
  return(list(opt.ret = ret, par=par.opt, value = ret$value))      
}

# A function similar to expand.grid, but different ordering of columns
make.grid.matrix = function(x=lapply(x.dim,function(n) 1:n),x.dim=NULL,n=NULL) {
restore.point("make.grid.matrix")
  #restore.point(""make.grid.matrix")
  
	if (!is.list(x)) {
  	# Simply a matrix
  	if (is.null(n) & is.null(x.dim)) {
    	return(x)
  	}
  	
		mat = matrix(NA,nrow=NROW(x)^n,ncol=n)
		for (i in 1:n) {
			mat[,i] = rep( rep(x,each=NROW(x)^(n-i)), times = NROW(x)^(i-1))
	  }
	  return (mat)
  } else {
	  n = length(x)
	  if (is.null(x.dim)) {
	  	x.dim = sapply(x,length)
  	}
	  mat = matrix(NA,nrow=prod(x.dim),ncol=n)
	  x.dim = c(1,x.dim,1,1)
		for (i in 1:n) {
			mat[,i] = rep( rep(x[[i]],each=prod(x.dim[(i+2):(n+2)])), times = prod(x.dim[1:i]))
	  }
	  return (mat)
	}
}



# Gives the corresponding rows for a permutated grid.matrix given a permutation
# x.perm of the elements of the original list x 
grid.matrix.permutation = function(x,perm.col) {
restore.point("grid.matrix.permutation")
  #restore.point(""grid.matrix.permutation")
  
	
	stopifnot(is.list(x))
	x.dim = sapply(x,length)

	x.dim.perm = x.dim[perm.col]
	gm = make.grid.matrix(x.dim=x.dim)
	rows.perm = rep(0,NROW(gm))
	
	nc = length(x.dim)
	if (nc==1) {
		return(1:NROW(gm))
	}
	for (mycol in 1:(nc-1)) {
		rows.perm = rows.perm + (gm[,perm.col[mycol]]-1)*prod(x.dim.perm[(mycol+1):nc])
	}
	rows.perm = rows.perm + gm[,perm.col[nc]]
	return(rows.perm)
}

# Example
 # x= list(c("A","B"),c("a","b"),c("0","1"))
 # perm.col = c(3,1,2)
 # x.perm =x[perm.col]
 
 # gm.x = make.grid.matrix(x)
 # gm.perm = make.grid.matrix(x.perm)
 
 # rows.perm = grid.matrix.permutation(x,perm.col)
 
 # cbind(gm.x,gm.perm[rows.perm,])



to.list = function(ob,just.return.if.list=TRUE,return.null=TRUE) {
  if (is.list(ob) & just.return.if.list)
    return(ob)
  li = list()
  li[[1]] = ob
  li
}
    
# My wrapper to the lattice function levelplot. Allows for some own color schemes
# The parameter focus specifies at which z range stronger color changes shall appear
sk.levelplot = function(x=NULL,y=NULL,z=NULL, xnames = NULL, ynames=NULL, grid.xyz = NULL,
                        col.scheme = "darkredgreen", na.col = NULL, at = NULL, at.scheme = "interval",
                        focus = 0,  cuts=15,col.regions=NULL, xlab=NULL,ylab=NULL, 
                        panel = panel.levelplot, zlim=NULL, reverse.colors=FALSE, ...) {

  #restore.point(""sk.levelplot")
  
  require(lattice)
  if (!is.null(z)) {
    z.vec = as.vector(z)
  } else {
    z.vec = grid.xyz[,3]
  }

  # Make cutpoints
  if (is.null(zlim))
    zlim = range(z.vec[is.finite(z.vec)])
    
  if (is.null(at)) {
    if (at.scheme == "interval") {
      at.rel = seq(0,1,length=cuts)
      at.rel = at.rel ^ (abs(focus)^sign(-focus))
      at = zlim[1] + diff(zlim)*at.rel
    } else if (at.scheme == "pretty") {
      at = pretty(z.vec,cuts)
    }
  }
  
  # Select colors
  num.color = cuts+1
  if (col.scheme == "grey") {
    col.regions = grey(seq(0,1,length=num.color))
    if (is.null(na.col)) {
      na.col="darkblue"
      na.col=hsv(1/6,1/2,1/2)
    }
  } else if (col.scheme == "darkredgreen" | col.scheme == "default" ) {
    col.regions = c(rainbow(num.color,start=5/6, end = 1/3,v = (0.3+0.7*(1:num.color)/num.color))) # From Magenta into green
    col.regions = c(rainbow(num.color,start=5/6, end = 1/4,v = (0.3+0.7*(1:num.color)/num.color))) # From Magenta into green

    if (is.null(na.col))
      na.col="grey"
  } else if (col.scheme == "own") {
    col.regions = col.regions[1:num.color]
  } else {
    stop(paste("col.scheme", col.scheme, " not implemented."))
  }
  
  if (reverse.colors) {
    col.regions=rev(col.regions)
  }
  
  # Set NA ROWS below the lowest value and assign NA color
  if (sum(!is.finite(z.vec)) > 0) {
    at = c(zlim[1]-1,at)
    at = c(zlim[1]-3,at)
    col.regions = c(na.col,col.regions)
    if (!is.null(z)) {
      z[!is.finite(z)] = zlim[1]-2
    } else {
      grid.xyz[!is.finite(z.vec),3] = zlim[1]-2
    }
  }
                               
  if (!is.null(z)) {
    if (!is.null(xnames))
      rownames(z) = xnames
    if (!is.null(ynames))
      colnames(z) = ynames
      
    row.values = 1:NROW(z); col.values = 1:NCOL(z)
    if (is.numeric(x))
      row.values=x
    if (is.numeric(y))
      column.values=y
    pl = levelplot(x=z,row.values = row.values, column.values = column.values,
              at=at,col.regions=col.regions,xlab=xlab,ylab=ylab,panel=panel,...)
  } else {
    if (is.null(xlab)) 
      xlab = names(grid.xyz)[1]
    if (is.null(ylab)) 
      ylab = names(grid.xyz)[2]
      
    colnames(grid.xyz) = c("x","y","z")
    grid.xyz = as.data.frame(grid.xyz)
    pl = levelplot(z~x*y,data=grid.xyz, at=at,col.regions=col.regions,xlab=xlab,ylab=ylab,
                   panel=panel,...)
  }
  
  print(pl)
}

# Helper function to discretize a continous distribution.
# F.vec is a finite vector containing the value of the cdf at M different points.
# The function generates an M dimension vector of probabilities summing up to 1
# that discretize the distribution
discretize.given.F.vec = function(F.vec) {
  M = NROW(F.vec)
  phi = numeric(M)
  phi[1] = F.vec[1] + 0.5 * (F.vec[2]-F.vec[1])
  phi[M] = 1-F.vec[M] + 0.5 * (F.vec[M]-F.vec[M-1])
  if (NROW(F.vec) > 2) {
    ind.left = 1:(M-2)
    ind.mid = ind.left +1
    ind.right = ind.left +2
    phi[ind.mid] = 0.5*((F.vec[ind.mid]-F.vec[ind.left]) + (F.vec[ind.right]-F.vec[ind.mid]))
  }    
  return(phi)  
}

# Calculate numerically the expected value given a cdf
calc.mean.from.F.fun = function(F.fun,x.min=0,x.max=Inf,abs.tol = 10^(-10),
                            x.seq=NULL,use.num.integrate=TRUE,...) {
  if (x.min >= 0 & use.num.integrate) {
    H.fun = function(x,...) {1-F.fun(x,...)} 
    mx = integrate(H.fun,lower=x.min,upper=x.max,abs.tol=abs.tol,...)$value
    return(mx)
  } else {
    if (is.null(x.seq)) {
      x.seq = seq(x.min,x.max,length=1000)
    }
    F.vec = F.fun(x.seq,...)
  	prob = discretize.given.F.vec(F.vec)
    mx = sum(prob*x.seq)
  }
  mx
} 


all.a.get.Ue.L = function(m,L) {
  if (is.null(m$LY)) {
    warning("all.a.get.Ue.L: m$LY not yet initialized. You have to call solve.all.a before.")
    return(NULL)
  }
  
  Ue.L.fun = function(a) {
    mat = m$LY[["e"]][[a]]
    if (!is.matrix(mat)) 
      return(NA)
    if (NROW(mat)==1) {
      if (L < mat[1,"L"]) {
        return(NA)
      } else {
        return(mat[1,"Ue"])
      }
    }
    approx(mat[,"L"],mat[,"Ue"],xout=L, method="linear",
           yleft=NA, yright=max(mat[,"Ue"]), rule = 1, f = 0, ties = "ordered")$y    
  }
  ret = sapply(1:m$nA, Ue.L.fun)
  ret
}
