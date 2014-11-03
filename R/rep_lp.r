
make.ac.for.lp = function(m,a,lp=m$lp, new.basis = TRUE, make.row.names=FALSE,
                          audit.prob = m$audit.prob, lambda = m$lambda) {
	
	#restore.point("make.ac.for.lp")	
	n = m$n; ny = m$ny; a.dim = m$a.dim;
	
	nc = n*ny+1;
	# Number of constraint rows for each type of constraint and indizes of those rows
	nr.ac =  sum(m$a.dim)- m$n;
	rows.ac = (1:nr.ac)+1;
		
	# Number of total non-zero entries for AC
	ne.ac = nr.ac * ny;
	ir = ic = v = rep(NA,ne.ac)
					
	# Action IC for every player i and all possible action profiles
	# For every signal we calculate m$phi.mat[,a]-m$phi.mat[,ak]

	if (make.row.names) {
		rowlab = rep("",nr.ac)
	}
	e0 = 0
	row = 2
	for (i in 1:n) {
		replies.i = get.replies.ind(m,a.ind=a,i=i)
		A.act = setdiff(replies.i,a)
		rowstart = row
		for (a.dev in A.act) {
			e.ind = e0+1:ny
			ir[e.ind]= row;
			ic[e.ind]= (1:ny)+(i-1)*ny+1;
			v[e.ind] = (1-audit.prob) * m$phi.mat[,a.dev]-m$phi.mat[,a]
			
			# Set RHS
			glp_set_row_bnds(lp,row,GLP_LO,m$g[a.dev,i]-m$g[a,i],Inf)

			row = row+1
			e0 = e0+ny
		}
		if (make.row.names) {
			rowlab[rowstart:(row-1)] = paste("AC i",i," ",m$lab.a[A.act],sep="")
		}
	}
	
	
	
	if (audit.prob > 0) {
  	# Assume all liquidity is given to a single player
  	i = which(lambda==1)
  	if (i == 1) {
    	rows =  2:(m$a.dim[1]-1)
  	} else {
    	row.start = 2 + sum(m$a.dim[1:(i-1)])-(i-1)
    	rows = row.start:(row.start+m$a.dim[i]-2)
  	}
  	ir.audit = rows
  	ic.audit = rep(1,NROW(ir.audit))
  	# The action constraint of player i will be relaxed by lambda.i*L
  	v.audit = rep(audit.prob,NROW(ir.audit))
	} else {
  	ir.audit = ic.audit = v.audit = NULL
	}
  	
	
	# Also create coefficients for Bmax (since those also depend on the action profile a
	# It repeats phi.mat[,a] for every player
	e.ind = 1:(nc-1)
	ir.B = rep(1,nc-1); ic.B = 2:nc; v.B = rep(m$phi.mat[,a], times=m$n)

	ir = c(ir.B,ir,ir.audit,m$lp.info$ir);
	ic=c(ic.B,ic,ic.audit,m$lp.info$ic);
	v=c(v.B,v,v.audit,m$lp.info$v);
	use.ir = (v!=0);

	# Load action constraints and the constraints for Bmax	
	glp_load_matrix(lp, NROW(ir[use.ir]), ir[use.ir], ic[use.ir], v[use.ir]);	

	# Set names
	if (make.row.names) {
		for (i in rows.ac) {
			glp_set_row_name(lp, i, rowlab[i]);
		}
	}
	#glp_print_prob(lp, "lp_prob.txt")
	if (new.basis) {
		glp_std_basis(lp);
	}
	
	return(lp)
}


# Modifies budget constraints such that at most a total amount B
# is allowed to be burned after any signal y
# Useful to look at renegotiation-proof equilibria or for characterizations
# of games with side payments but without money burning
set.money.burning.constraint = function(m,B = Inf, lp = m$lp) {
	
	#restore.point("make.default.lp.for.m")		

	require(glpkAPI)
	nc = n*ny+1;
	nr.pc = m$ny; nr.ac =  sum(m$a.dim)- m$n; nr.bc = m$ny;
	rows.bc = ((nr.ac+nr.pc+1):(nr.ac+nr.pc+nr.bc))+1;
		
	# Budget constraints satisfy  0 <= sum(p.i(y)) <= B
	if (is.finite(B)) {
  	for (r in rows.bc) {
  		glp_set_row_bnds(lp,r,GLP_DB,0,B)
  	}
	} else {
  	for (r in rows.bc) {
  		glp_set_row_bnds(lp,r,GLP_LO,0,Inf)
  	}
  }
}


# The linear problem is more compact if all liquidity is given just to one player
# i denotes the player who gets the liquidity

# Attention! Do not use i as some counter in the code
make.default.lp.for.m = function(m, i=m$n, audit.prob = m$audit.prob) {
	
	#restore.point("make.default.lp.for.m")		

	library(glpkAPI)
	n = m$n; ny = m$ny; a.dim = m$a.dim;
	
	nc = n*ny+1;
	
	# Number of constraint rows for each type of constraint and indizes of those rows
	# We only add explicit payment constraints for player i. The other payment constraints will be 
	# modelled by a zero upper bound on p.j(y), i.e. players can only receive money from player i
	
	#nr.pc = ny*m$n; nr.ac =  sum(m$a.dim)- m$n; nr.bc = m$ny;
	nr.pc = ny; nr.ac =  sum(m$a.dim)- m$n; nr.bc = m$ny;

	nr= nr.ac+nr.pc+nr.bc +1;
	rows.ac = (1:nr.ac)+1;
	rows.pc = ((nr.ac+1):(nr.ac+nr.pc))+1;
	rows.bc = ((nr.ac+nr.pc+1):(nr.ac+nr.pc+nr.bc))+1;
		
	# Number of total non-zero entries
	ne.pc = nr.pc * 2;
	#ne.ac = nr.ac * ny;
	ne.bc = nr.bc * 2;
	
	ne = ne.pc + ne.bc
	
	# Indices for sparse matrix
	ir = ic = v = rep(NA,ne) 

	
	# PC for player i and every signal
	# The first column (L) of the payment constraints are 1
	# Afterwards we have an ny x ny idenity matrix starting at the place of player i

	e0 = 1
	e.ind = e0:(e0+nr.pc-1)
	# Coefficient before L in PC
	ir[e.ind]=rows.pc; ic[e.ind]=1; v[e.ind] = -1
	# Coefficient before L in PC
	e.ind = e.ind+nr.pc
	#ir[e.ind]=rows.pc; ic[e.ind]=(1:nr.pc)+1; v[e.ind] = 1
	ir[e.ind]=rows.pc; ic[e.ind]=(1:nr.pc)+1 + (ny*(i-1)); v[e.ind] = 1

		
# 	e0 = 1
# 	e.ind = e0:(e0+nr.pc-1)
# 	# Coefficient before L in PC
# 	ir[e.ind]=rows.pc; ic[e.ind]=1; v[e.ind] = -1
# 	# Coefficient before L in PC
# 	e.ind = e.ind+nr.pc
# 	#ir[e.ind]=rows.pc; ic[e.ind]=(1:nr.pc)+1; v[e.ind] = 1
# 	ir[e.ind]=rows.pc; ic[e.ind]=(1:nr.pc)+1 + (ny*(i-1)); v[e.ind] = 1
	
	# BC for every signal y
	# Basically n identity matrixes of size ny*ny next to each other
	# The first column for L is 0 everywhere
	e0 = ne.pc+1;
	e.ind = e0:(e0+ny-1)
	for (j in 1:n) {
		ir[e.ind]=rows.bc; ic[e.ind]=(1:ny)+(j-1)*ny+1; v[e.ind] = 1
		e.ind = e.ind+ny
	}
	
	# RHS and DIR
	

	# Dimnames
	#lab.pc = paste("PC ", rep(1:n,each=ny),"(",rep(m$lab.y[1:ny],times=n),")",sep="")
	lab.ac  = NULL
	for (j in 1:m$n) {
		lab.ac = c(lab.ac,paste("AC ",j,".",1:(m$a.dim[j]-1),sep=""))
	}
	lab.pc = paste("PC ",m$lab.y[1:ny],sep="")
	rowlab = c("Bmax",lab.ac, lab.pc, paste("BC ",m$lab.y[1:ny],sep=""))
	str = NULL
	for (j in 1:m$n) str = c(str,paste("p",j,"_",m$lab.y,sep=""));
	collab = c("L",str)
	
	#mat = simple_triplet_matrix(ir, ic, v, nrow = nr, ncol = nc, dimnames = list(rowlab,collab))
	
	#################################################################################################
	# Initialize lp using GLPK
	#################################################################################################		
	
	lp <- create.lp.prob()
		
	glp_add_cols(lp, nc)
	glp_add_rows(lp, nr) 

	#glp_load_matrix(lp, NROW(ir), ir, ic, v);
		
	glp_set_prob_name(lp, "LP Stub")
	glp_set_obj_dir(lp, GLP_MIN)
		
	# Set default bounds for variables 
	# GLP_FR    -8 <  x  <  +8    Free (unbounded) variable
	# GLP_LO    lb <= x  <  +8    Variable with lower bound
	# GLP_UP    -8 <  x  <= ub    Variable with upper bound
	# GLP_DB    lb <= x  <= ub    Double-bounded variable
	# GLP_FX    lb =  x  =  ub    Fixed variable
	
	# For L: Set lower bound 0 and no upper bound
	glp_set_col_bnds(lp,1,GLP_LO,0,Inf)

	# Set upper bound of 0 for p.j(y) for j not i
	# 
  for (j in setdiff(1:m$n,i)) {
	  for (y in 1:ny) {
	  	glp_set_col_bnds(lp,1+(j-1)*ny+y,GLP_UP,-Inf,0)
  	}
	}  
	# Set p.i(y) unbounded
	for (y in 1:ny) {
		glp_set_col_bnds(lp,1+(i-1)*ny+y,GLP_FR,-Inf,Inf)
  }

	# Set row bounds and directions 
# 	btype[dir=="<="] = GLP_UP
#  	btype[dir==">="] = GLP_LO
#  	btype[dir=="free"] = GLP_FR

	
	# Bmax is initially free
	glp_set_row_bnds(lp,1,GLP_FR,0,0)
	
	# Action constraint bounds must be set latter
	# glp_set_row_bnds(lp,i,GLP_LO,lp,m$g[a,i]-m$g[ak,i],Inf)

	# Payment Constraints satisfy <= 0 
	for (r in rows.pc) {
		glp_set_row_bnds(lp,r,GLP_UP,-Inf,0)
	}
	
	
	# Budget constraints satisfy >= 0 
	for (r in rows.bc) {
		glp_set_row_bnds(lp,r,GLP_LO,0,Inf)
	}

	# Set column and row names	
	for (j in 1:nc) {
		glp_set_col_name(lp, j, collab[j]);
	}
	for (r in 1:nr) {
		glp_set_row_name(lp, r, rowlab[r]);
	}
	
	#glp_print_prob(lp, "lp_prob.txt")	

	lambda = rep(0,m$n)
	lambda[i] = 1
	lp.info = list(ir=ir,ic=ic,v=v,nr=nr,nc=nc, a=NULL, lambda=lambda, i = i)
			
	return(list(lp=lp,lp.info=lp.info))
}



make.sym.default.lp.for.m = function(m) {
	
	#assign.list(VAR.STORE[["make.default.lp.for.m"]]$var)		

	library(glpkAPI)
	n = m$n; ny = m$ny; a.dim = m$a.dim;
	
	nc = n*ny+1;
	# Number of constraint rows for each type of constraint and indizes of those rows
	nr.pc = n*ny; nr.ac =  sum(m$a.dim)- m$n; nr.bc = m$ny;
	nr= nr.ac+nr.pc+nr.bc +1;
	rows.ac = (1:nr.ac)+1;
	rows.pc = ((nr.ac+1):(nr.ac+nr.pc))+1;
	rows.bc = ((nr.ac+nr.pc+1):(nr.ac+nr.pc+nr.bc))+1;
		
	# Number of total non-zero entries
	ne.pc = nr.pc * 2;
	#ne.ac = nr.ac * ny;
	ne.bc = nr.bc * 2;
	
	ne = ne.pc + ne.bc
	
	# Indices for sparse matrix
	ir = ic = v = rep(NA,ne) 

	
	# PC for every player and every signal
	# The first column (L) of the payment constraints are 1/n
	# Afterwards we have an idenity matrix
	
	e0 = 1
	e.ind = e0:(e0+nr.pc-1)
	# Coefficient before L in PC
	ir[e.ind]=rows.pc; ic[e.ind]=1; v[e.ind] = -1/m$n
	# Coefficient before L in PC
	e.ind = e.ind+nr.pc
	ir[e.ind]=rows.pc; ic[e.ind]=(1:nr.pc)+1; v[e.ind] = 1
	
	# BC for every signal y
	# Basically n identity matrixes of size ny*ny next to each other
	# The first column for L is 0 everywhere
	e0 = ne.pc+1;
	e.ind = e0:(e0+ny-1)
	for (i in 1:n) {
		ir[e.ind]=rows.bc; ic[e.ind]=(1:ny)+(i-1)*ny+1; v[e.ind] = 1
		e.ind = e.ind+ny
	}
	
	# RHS and DIR
	

	# Dimnames
	lab.pc = paste("PC ", rep(1:n,each=ny),"(",rep(m$lab.y[1:ny],times=n),")",sep="")
	rowlab = c("Bmax",rep("AC",nr.ac), lab.pc, paste("BC ",m$lab.y[1:ny],sep=""))
	str = NULL
	for (i in 1:m$n) str = c(str,paste("p",i,"_",m$lab.y,sep=""));
	collab = c("L",str)
	
	#mat = simple_triplet_matrix(ir, ic, v, nrow = nr, ncol = nc, dimnames = list(rowlab,collab))
	
	#################################################################################################
	# Initialize lp using GLPK
	#################################################################################################		
	
	lp <- create.lp.prob()
		
	glp_add_cols(lp, nc)
	glp_add_rows(lp, nr) 

	#glp_load_matrix(lp, NROW(ir), ir, ic, v);
		
	glp_set_prob_name(lp, "LP Stub")
	glp_set_obj_dir(lp, GLP_MIN)
		
	# Set default bounds for variables 
	# GLP_FR    -8 <  x  <  +8    Free (unbounded) variable
	# GLP_LO    lb <= x  <  +8    Variable with lower bound
	# GLP_UP    -8 <  x  <= ub    Variable with upper bound
	# GLP_DB    lb <= x  <= ub    Double-bounded variable
	# GLP_FX    lb =  x  =  ub    Fixed variable
	
	# For L: Set lower bound 0 and no upper bound
	glp_set_col_bnds(lp,1,GLP_LO,0,Inf)

	# For p.i(y) set unbounded
  for (j in 2:nc) {
	  glp_set_col_bnds(lp,j,GLP_FR,-Inf,Inf)
	}  

	# Set row bounds and directions 
# 	btype[dir=="<="] = GLP_UP
#  	btype[dir==">="] = GLP_LO
#  	btype[dir=="free"] = GLP_FR

	
	# Bmax is initially free
	glp_set_row_bnds(lp,1,GLP_FR,0,0)
	
	# Action constraint bounds must be set latter
	# glp_set_row_bnds(lp,i,GLP_LO,lp,m$g[a,i]-m$g[ak,i],Inf)

	# Payment Constraints satisfy <= 0 
	for (i in rows.pc) {
		glp_set_row_bnds(lp,i,GLP_UP,-Inf,0)
	}
	# Budget constraints satisfy >= 0 
	for (i in rows.bc) {
		glp_set_row_bnds(lp,i,GLP_LO,0,Inf)
	}

	# Set column and row names	
	for (j in 1:nc) {
		glp_set_col_name(lp, j, collab[j]);
	}
	for (i in 1:nr) {
		glp_set_row_name(lp, i, rowlab[i]);
	}
	
	#glp_print_prob(lp, "lp_prob.txt")	

	lp.info = list(ir=ir,ic=ic,v=v,nr=nr,nc=nc,rowlab=rowlab,collab=collab, a=NULL)

			
	return(list(lp=lp,lp.info=lp.info))
}

