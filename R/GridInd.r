# Consider a cartesian grid X = X1 x ... x Xm where Xk is some list.
# For example, X could be a set of action profiles and Xk is the set of actions of player i
# We sometimes may want to assign to every element of X one or more values, e.g. payoffs for every player

# It is then sometimes fast and convenient to store these values in a matrix where every row corresponds to
# one element of X, i.e. the matrix has n.X = n.X1 * ... * n.Xm rows.

# This package provides some useful index functions for such vector grids


#' Generates a GridIndex
GridInd = function(val=NULL,m=NULL,dim=NULL, names=NULL) {
  gi = list()
  if (!is.null(m)) {
    gi$val = list()
    for (i in 1:m) {
      gi$val[[i]] = val
    }
  } else if (!is.null(dim) & is.null(val)) {
    m = length(dim)
    gi$val = list()
    for (i in 1:m) {
      gi$val[[i]] = 1:dim[m]
    }
  } else {
    gi$val = val
  }
  
  gi$m   = length(gi$val)
  gi$dim = sapply(gi$val,length)
  gi$n   = prod(gi$dim)
  
  gi$shift   = rev(cumprod(rev(gi$dim))) / gi$dim
	gi$modulus = c(prod(gi$dim),gi$shift[-gi$m])
	
	gi$names = names
  class(gi) = c("GridInd","list")

	return(gi)
}


examples.GridInd = function() {
 gi = GridInd(1:2,m=3)
 gi
 grid.matrix(gi)
 v.ind.to.x.ind(gi)

 make.replies.gi(gi,k.i=c(1))
 cgi = make.child.gi(gi,k=c(2,3))
 grid.matrix(cgi)
 cbind(cgi$cp.ind,grid.matrix(gi))

 rep.gi = make.replies.gi(gi,k.i=c(1,2))
 rep.gi
 cbind(rep.gi$cp.ind,grid.matrix(gi))
}


# x value to an index

x.to.x.ind = function(gi,x) {
  if (is.matrix(x)) {
    x.ind = matrix(n,nrow=NROW(x),ncol=NCOL(x))
    for (k in 1:gi$m)
      x.ind[,k] = match(x[,k],gi$val[[k]])
    return(x.ind)
  }
  x.ind = numeric(length(x))
  for (k in 1:gi$m)
    x.ind[k] = match(x[k],gi$val[[k]])
  return(x.ind)
}



# matrix of x.indices to x.values
x.ind.to.x = function(gi,x.ind,k=1:gi$m) {
  if (is.matrix(x.ind)) {
    x = matrix(NA,nrow=NROW(x),ncol=NCOL(x))
    for (k.act in k)
      x[,k.act] = gi$val[[k.act]][x.ind]
  } else if (length(k)==1) {
    x = gi$val[[k]][x.ind]
  } else {
    x = numeric(length(x))
    for (i in 1:length(x)) {
      k.act = k[i]
      x[i] = gi$val[[k.act]][x.ind[i]]
    }
  }
  return(x.ind)
}


# x value to an index
x.to.x.ind = function(gi,x,k=1:gi$m) {
  if (is.matrix(x)) {
    x.ind = matrix(NA,nrow=NROW(x),ncol=NCOL(x))
    for (k.act in k)
      x.ind[,k.act] = match(x[,k.act],gi$val[[k.act]])
  } else if (length(k)==1) {
    x.ind = match(x,gi$val[[k]])
  } else {
    x.ind = numeric(length(x))
    for (i in 1:length(x)) {
      k.act = k[i]
      x.ind[i] = match(x[i],gi$val[[k.act]])
    }
  }
  return(x.ind)
}

x.to.v.ind = function(gi,x) {
  x.ind = x.to.x.ind(gi,x)
  x.ind.to.v.ind(gi,x.ind)
}

x.ind.to.v.ind = function(gi,x.ind) {
  if (is.matrix(x.ind)) {
  	v = rep(1,NROW(x.ind))
  	for (k in 1:gi$m) {
  		v = v + (x.ind[,k]-1)*gi$shift[k]
  	}
  	return(v)
  } else {
    return(sum((x.ind-1)*gi$shift)+1)
  }
}

v.ind.to.x.ind = function(gi,v.ind=1:gi$n, k=1:gi$m,as.matrix = length(k)>1) {
  
  #restore.point("v.ind.to.x.ind")
  n.k = length(k)
  if (as.matrix) {
    x.ind = matrix(0,NROW(v.ind),n.k)
  	for (i in 1:n.k) {
  		x.ind[,i] = ceiling((((v.ind-1) %% gi$modulus[k[i]])+1) / gi$shift[k[i]])
  	}
  	return(x.ind)
  } else {
    if (n.k > 1) {
      x.ind = numeric(n.k)
  	  for (i in 1:n.k) {
    		x.ind[i] = ceiling((((v.ind-1) %% gi$modulus[k[i]])+1) / gi$shift[k[i]])
    	}
    	return(x.ind)
  	} else {
      x.ind = ceiling((((v.ind-1) %% gi$modulus[k])+1) / gi$shift[k])
    	return(x.ind)
    }      	
	}
}


"[.GridInd" <- function (vl,i,j,drop) {
  stop()
#   
#   ijcommas = length(sys.call()) -3 - (!missing(empty)) - (!missing(drop))
#   mycall = sys.call()
#   
#   #print(as.character(mycall))
#   #print(ijcommas)
#   if (missing(empty)) empty = NA
#   if (missing(drop)) drop = TRUE    
#    
#   if (missing(i) & missing(j))  {
#     # Called [,]
#     if (ijcommas>0)
#       return(as.matrix.VectorListInd(vl,vec=1:length(vl),empty=empty))
#     # Called []
#     return(1:length(vl))
#   }  
# 
#   if (missing(j)) {
#     if (ijcommas == 0 | (NROW(i) == 1 & drop==TRUE)) {
#       return(rows.to.v.ind(vl,rows=i))
#     } else {
#       return(as.matrix.VectorListInd(vl,vec = rows.to.v.ind(vl,rows=i), rows=i,empty=empty))
#     }
#   }
#   if (!missing(i) & !missing(j)) {
#     if (NROW(i) == 1) {
#       return(rows.to.v.ind(vl,rows=i)[j])
#     }
#   }
#   stop("Column indexing not yet fully implemented for VectorListInd")
}



v.ind.to.x = function(gi,v.ind=1:gi$n,k=1:gi$m, as.matrix = TRUE) {
  x.ind = v.ind.to.x.ind(gi,v.ind,k,as.matrix)
  x.ind.to.x(gi,x.ind,k)    
}

gi.x.ind.k = function(gi,k) {
  if (length(k)>1) 
    stop("only works for a single index k")
  x = 1:gi$dim[k]
  mydim = c(1,gi$dim,1,1)
  rep(rep(x,each=prod(mydim[(k+2):(gi$m+2)])), times = prod(mydim[1:k]))
}  


# Rough idea: A matrix of replies of player i for each a_i, i.e. action profile of
# other players
# Details: Not only generate the matrix but also a child grid index
# Assume that only a subset of the m vectors will vary
make.replies.gi = function(gi, k.i) {
  k_i = (1:gi$m)[-k.i]
  cgi = make.child.gi(gi,k_i)
  cgi$k.i = k.i
  cgi$reply.mat = matrix(order(cgi$cp.ind),nrow=cgi$n,byrow=TRUE)
  return(cgi)
}


# Assume we want an index in which every row corresponds to a profile a_k
# A child index is a grid index for this subgrid 
# In addition it has a vector cp.ind that gives the child row for every row of the parent index

make.child.gi = function(gi, k) {
  
  #restore.point("make.child.gi")

  cgi = GridInd(gi$val[k])
  # Generate index with cgi$parent$n rows. Each element gives the corresponding index number in cgi
  # This means vec[row] will be the index in cgi where the columns k are the same as in gi
  cgi$cp.ind = rep(1,gi$n)
  for (i in 1:length(k)) {
    x.ind = gi.x.ind.k(gi,k[i])
    cgi$cp.ind = cgi$cp.ind + (x.ind-1) * cgi$shift[i]
  }
  cgi$parent = gi

  return(cgi)
}



#Some simple tests
grid.matrix  = function(gi) {
   make.grid.matrix(x=gi$val)
 }
