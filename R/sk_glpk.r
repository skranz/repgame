
#glp_print_sol(lp, "lp_sol.txt")	


library(glpkAPI)

# Always use this function to create a new lp.prob
# the pointer is stored so that one can free the memory at any time

	

create.lp.prob = function(nr=0,nc=NULL) {
		lp <- glp_create_prob()
		#if (nc>0) glp_add_cols(lp, nc);
		#if (nr>0) glp_add_rows(lp, nr);

		#lp <- glp_create_prob()
		REP.CNT$LP_PROBS <<- c(REP.CNT$LP_PROBS,lp)
		lp
}

remove.LP_PROBS = function() {
	if (!is.null(REP.CNT$LP_PROBS)) {
		for (i in 1:length(REP.CNT$LP_PROBS)) {
		  glp_delete_prob(REP.CNT$LP_PROBS[[i]])
		}
	}
	REP.CNT$LP_PROBS <<- NULL
}
#remove.LP_PROBS()

remove.lp = function(lp) {
	try(glp_delete_prob(lp))
}

lp.get.col.shadow.price = function(lp,col) {
	glp_get_col_dual(lp,col)
}


lp.set.name = function(lp,lab) {
	#name.lp(lp, lab)
	glp_set_prob_name(lp,lab)
}


lp.set.obj.coef = function(lp,obj.coef) {
	for (j in 1:length(obj.coef)) {
	 glp_set_obj_coef(lp, j, obj.coef[j])
	}
}


make.glpk.lp = function(obj.coef=NULL, con=NULL, dir=NULL, rhs=NULL, obj.max = TRUE, 
		             bounds = NULL, row.lab=NULL, col.lab=NULL, prob.name=NULL, prob = NULL, lp=NULL,
		             set.names = TRUE
) {
  restore.point("make.glpk.lp")
  #restore.point("make.glpk.lp")
  
  
	library(glpkAPI)
	library(slam)
	
	if (!is.null(prob)) {
		obj.coef=prob$obj.coef; con=prob$con; dir= prob$dir;
		rhs = prob$rhs; obj.max = prob$obj.max; bounds = prob$bounds;
		row.lab = prob$row.lab;	col.lab = prob$col.lab; prob.name=prob.name
	}
	
	nc = NCOL(con)
	nr = NROW(con)
	
	if (is.null(lp)) {
		lp <- create.lp.prob()
		
		glp_add_cols(lp, nc)
		glp_add_rows(lp, nr) 
	}
	
	if (is.null(prob.name)) {
		glp_set_prob_name(lp, prob.name)
	}
	
	if (obj.max) {	
		glp_set_obj_dir(lp, GLP_MAX)
	} else {
		glp_set_obj_dir(lp, GLP_MIN)
	}
	
	if (!is.null(obj.coef)) {
		for (j in 1:nc) {
		 glp_set_obj_coef(lp, j, obj.coef[j])
		}
	}

	    
# GLP_FR    -8 <  x  <  +8    Free (unbounded) variable
# GLP_LO    lb <= x  <  +8    Variable with lower bound
# GLP_UP    -8 <  x  <= ub    Variable with upper bound
# GLP_DB    lb <= x  <= ub    Double-bounded variable
# GLP_FX    lb =  x  =  ub    Fixed variable
  
	if (is.null(bounds)) {
		lb = rep(0,nc)
		ub = rep(0,nc)
		btype = rep(GLP_FR,nc)
	}	else {
		lb = bounds$lower
 		ub = bounds$upper 		 		
		
    #lb = bounds$lower$val
		#ub = bounds$upper$val 		 		
		
    btype = rep(GLP_DB,nc)
 		btype[lb == ub ] = GLP_FX
 		btype[!is.finite(lb) & !is.finite(ub)] = GLP_FR
 		btype[is.finite(lb) & !is.finite(ub)] = GLP_LO
 		btype[!is.finite(lb) & is.finite(ub)] = GLP_UP
	}
	
	#Set bounds of the variables
  #glp_set_col_bnds(lp, j, type, lb, ub)
  #mapply(glp_set_col_bnds, j = 1:nc, type=btype, lb=lb,ub=ub, MoreArgs=list(lp=lp))
  for (j in 1:nc) {
	  glp_set_col_bnds(lp,j,btype[j],lb[j],ub[j])
	}  

	
	lb = ub = rep(0,nr)
	btype = rep(GLP_FX,nr)
 	btype[dir=="<="] = GLP_UP
 	btype[dir==">="] = GLP_LO
 	btype[dir=="free"] = GLP_FR
	ub[dir=="<="] = rhs[dir=="<="]
	lb[dir==">="] = rhs[dir==">="]
	
	for (i in 1:nr) {
	  glp_set_row_bnds(lp,i,btype[i],lb[i],ub[i])
 	}
	
 	con.tri <- as.simple_triplet_matrix(con)
	glp_load_matrix(lp, length(con.tri$i), con.tri$i, con.tri$j, con.tri$v);
	
	
	if (set.names) {
		if (!is.null(col.lab)) {
			for (j in 1:nc) {
				glp_set_col_name(lp, j, col.lab[j]);
			}
		}
		
		if (!is.null(row.lab)) {
			for (i in 1:nr) {
				glp_set_row_name(lp, i, row.lab[i]);
			}
		}
	}
		
	return(lp)
}

solve.glpk.lp = function(lp, sink_file = REP.CNT$SINK_FILE, retry.with.standard.basis = TRUE, 
         warning.if.no.solution = TRUE, retry.with.dual = retry.with.standard.basis,
         call.lab="", should.always.solve = TRUE, max.iteration = REP.CNT$LP.MAX.ITERATION) {
  
  restore.point("solve.glpk.lp")
  #restore.point("solve.glpk.lp")         
           
	LAST_SOLVE_LP <<- lp
	nr = glp_get_num_rows(lp)
	nc = glp_get_num_cols(lp)
	
	setSimplexParmGLPK(IT_LIM,max.iteration)
	#glp_set_int_parm(lp, glp_K_ITLIM, max.iteration)
	ret = list();
	print(paste("solve.lp ", call.lab))
	if (!is.null(sink_file)) {
		sink(sink_file)
		glp_simplex(lp);
		sink(NULL)
	}	else {
		glp_simplex(lp);
		#flush.console()
	}
	ret$glp.status = 	glp_get_status(lp)
  ret$status.lab = status_codeGLPK(ret$glp.status)
	
  ret$status = (ret$glp.status != GLP_OPT)*1
  
  worked = NULL
  worked = if (ret$status == 0) worked = "direct"
  
  if (ret$status != 0 & retry.with.standard.basis) {
  	REP.CNT$LP.RETRY.COUNT <<- REP.CNT$LP.RETRY.COUNT+1
  	
		setSimplexParmGLPK(IT_LIM,max.iteration)
  	#glp_set_int_parm(lp, glp_K_ITLIM, max.iteration) 
	  glp_std_basis(lp);
	  glp_simplex(lp);
		#flush.console()

	  ret$status = 	glp_get_status(lp)
  	ret$status = (ret$glp.status != GLP_OPT)*1
  	if (ret$status == 0) worked = "basis"
	}
	
  if (ret$status != 0 & retry.with.dual) {
    
  	REP.CNT$LP.RETRY.COUNT <<- REP.CNT$LP.RETRY.COUNT+1
  	
  	setSimplexParmGLPK(IT_LIM,max.iteration)
		#glp_set_int_parm(lp, glp_K_ITLIM, max.iteration)  	
  	setSimplexParmGLPK("meth",GLP_DUAL)
		#glp_set_int_parm(lp, glp_K_DUAL, 1)  	  
	  
		
		glp_std_basis(lp);
	  glp_simplex(lp);
		#flush.console()
	  
		setSimplexParmGLPK("meth",GLP_DUAL)
  	#glp_set_int_parm(lp, glp_K_DUAL, 0)   

	  ret$status = 	glp_get_status(lp)
  	ret$status = (ret$glp.status != GLP_OPT)*1
  	if (ret$status == 0) worked = "dual"
	}
	
	ret$glp.status = 	glp_get_status(lp)
  ret$status.lab = status_codeGLPK(ret$glp.status)	
	if (ret$status == 0) {
	  if (worked != "direct") {
  	  if (worked == "basis")
  		  wstr = paste("solve.glpk.lp from ", call.lab,
			          ": had to reinitilaize basis to solve the problem.")
			          
  	  if (worked == "dual")
  		  wstr = paste("solve.glpk.lp from ", call.lab,
			          ": had to use dual simplex to solve the problem.")
			myglpk.do.lp.warning(ret=ret,lp=lp, wstr=wstr, warning.if.no.solution=TRUE, should.always.solve=FALSE)
    }
  }	else {
  	wstr = paste("solve.glpk.lp from ", call.lab,
			          ": could not solve the problem.")
		myglpk.do.lp.warning(ret=ret,lp=lp, wstr=wstr, warning.if.no.solution, should.always.solve)	
  }			
  
	ret$solution     = mapply(glp_get_col_prim,j=1:nc,MoreArgs=list(lp=lp))
	ret$shadow.price = mapply(glp_get_row_dual,i=1:nr,MoreArgs=list(lp=lp))
	ret$objval       = glp_get_obj_val(lp)

	ret
}

myglpk.do.lp.warning = function(ret,lp, wstr, warning.if.no.solution = TRUE, should.always.solve=TRUE) {
  if (warning.if.no.solution | (ret$glp.stat != 183 & ret$glp.stat != 185) ) {
    warning(wstr)
    print("")
    #print(paste("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
    print(paste("warning solve.lp ",status_codeGLPK(ret$glp.status), wstr))
    print("")
  	SOLVE.LP.WARNINGS = c(REP.CNT$SOLVE.LP.WARNINGS, wstr)
  	fname = paste("lp_sol_warn",NROW(REP.CNT$SOLVE.LP.WARNINGS),".txt",sep="")
  	glp_print_sol(lp, fname)
  	flush.console()
	}
	if (should.always.solve) {
	  ERROR.GLPK.SOLVE <<- TRUE			  
	}		            
}


info.glpk.lp = function(lp, fname="lp_sol.txt") {
	glp_print_sol(lp, fname)
}	

bounds.to.Rglpk.bounds = function(bounds,nc = NULL)	{
		if (bounds==NULL) {
			return(NULL)
		} 
		if (is.null(nc)) {
			nc =length(bounds(lower))
		}
		bounds <- list(lower = list(ind = c(1:nc), val = bounds$lower),
    	             upper = list(ind = c(1:nc), val = bounds$upper))
    bounds
}

glpk.set.max.time = function(lp,max.time) {
   #glp_K_ITLIM, glp_K_ITCNT:
   #glp_set_int_parm(lp, glp_K_TMLIM, max.time)
	 setSimplexParmGLPK(TM_LIM,max.time)
}	
