# Functions for plotting

#' Plots payoffs as function of the discount factor
get.payoff.at.delta = function(m,delta= seq(0,1, length=101), payoff.var = "Ue") {
  restore.point("get.payoff.at.delta")
  # restore.point("get.payoff.at.delta")

  mat = m$opt.mat
  if (m$sol.type=="imp") {
    ret = approx(x = mat[,"delta"], y = mat[,payoff.var], xout=delta, method="linear",
               yleft=Inf, yright=mat[NROW(mat),payoff.var], rule = 1, f = 0, ties = "ordered")

  } else {
    mat = m$opt.mat
    ret = approx(x = mat[,"delta"], y = mat[,payoff.var], xout=delta, method="constant",
               yleft=Inf, yright=mat[NROW(mat),payoff.var], rule = 1, f = 0, ties = "ordered")
  }
  return(ret$y)
}
get.payoffs.from.model.list = Vectorize(get.payoff.at.delta, vectorize.args="m",SIMPLIFY=TRUE)

#' Creates a levelplot of payoffs from a model.list with parameters or the discount factor on the x and y axes
levelplot.payoff.compstat = function(m.list, par ,  xvar = "par1", yvar = "delta", payoff.var="Ue", delta = NULL,
   par.labels= round(par,2),identify = NULL, col.scheme = "grey", main = NULL,use.rows = NULL,...) {
	
  library(lattice)
  
	#restore.point("levelplot.payoff.compstat")
	
  if (!is.matrix(par)) {
    par = matrix(par,NROW(par),1)
    if (xvar != "delta") {
      colnames(par) = xvar
    } else {
      colnames(par) = yvar
    }
  }

  
  transpose = FALSE
  delta.is.par = FALSE
  if (is.null(delta)) {
    if (xvar == "delta" | yvar == "delta") {
	    delta = seq(0,1,by=0.001)
    } else {
      delta = 0
    }
  }
  if (!(xvar == "delta" | yvar == "delta")) {
    if (NROW(delta)>1) {
      warning("You gave a vector of delta, but neither xvar nor yvar is 'delta'. I pick first element of delta.")
      delta = delta[1]
    }
  }

  payoffs = get.payoffs.from.model.list(m.list,delta=delta,payoff.var=payoff.var)

  if (is.null(main)) {
    if (NROW(delta) > 1) {
      main = paste(payoff.var)
    } else {
      main = paste(payoff.var, " delta = ", round(delta,4))
    }
  }          
     
  if (xvar == "delta" | yvar == "delta") {
    par.col = xvar
    if (xvar == "delta") par.col = yvar
    grid.xy = make.grid.matrix(x=list(delta=delta, par.col=par[,par.col]))
    payoffs = as.vector(t(payoffs))
    if (!is.null(use.rows)) {
      grid.xy = grid.xy[use.rows,,drop=FALSE]
      payoffs = payoffs[use.rows,,drop=FALSE]
    }
    if (xvar == "delta") {
      grid.xyz = data.frame(grid.xy[,1],grid.xy[,2],payoffs)
    } else {
      grid.xyz = data.frame(grid.xy[,2],grid.xy[,1],payoffs)
    }
    colnames(grid.xyz)=c(xvar,yvar,payoff.var)
    sk.levelplot(grid.xyz = grid.xyz, main = main,
             panel = panel.levelplot,col.scheme = col.scheme,...)
  } else {
    if (!is.null(use.rows)) {
      par = grid.xy[use.rows,,drop=FALSE]
      payoffs = payoffs[use.rows,,drop=FALSE]
    }
    grid.xyz = data.frame(par[,xvar],par[,yvar],payoffs)
    colnames(grid.xyz)=c(xvar,yvar,payoff.var)
    sk.levelplot(grid.xyz = grid.xyz, main = main,
             panel = panel.levelplot,col.scheme = col.scheme,...)
  }

	
  if (isTRUE(identify==FALSE))
    identify = NULL
  
  if (!is.null(identify)) {
    lab = rep("",NROW(grid.xyz))
    trellis.focus("panel", 1, 1)
    
    a.vec = NULL
    dev.new = NULL
        
    while(TRUE) {
      ind = panel.identify(n=1,labels = lab)
      if (NROW(ind)<1) break;
      pr = grid.xyz[ind,,drop=FALSE]
      lpoints(grid.xyz[ind,1],grid.xyz[ind,2],col="blue",pch=3)
      colnames(pr) = c(xvar,yvar,payoff.var)
      print(pr)         
      flush.console()
    }
  }
  
  #invisible(return(grid.xyz))   
}


#' Shows a levelplot that illustrates the solution of a simple two player repeated game with imperfect monitoring
levelplot.rep = function(m, z=m$G, identify = "Ue",col.scheme = "grey", main = "",iter.lim = 20000,na.col=NULL,...) {
	
	#restore.point("levelplot.rep")
  
	if (m$sol.type != "imp") {
  	warning("levelplot.rep so far only implemented for games with imperfect public monitoring")
  	return(NULL)
	}
	#sk.levelplot(grid.xyz = data.frame(m$action.val.mat, z), m=m, main = main,
  #           panel = panel.levelplot.repgame,col.scheme = col.scheme)

	sk.levelplot(grid.xyz = data.frame(m$action.val.mat, z), m=m, main = main,
             panel = panel.levelplot.repgame,col.scheme = col.scheme,na.col=na.col,...)

  
  col.a = rainbow(100)
  col.a = col.a[
    as.vector(sapply(c(1,11,5,15,3,8,13,18,2,7,12,17,6,19,10,20),function(offset) 0:4*20+offset))
  ]
  if (is.null(m$LUe)) {
    m$LUe = as.list(rep(NA,m$nA))
  }
  xli <<- yli <<- c(Inf,-Inf)
  xli[1] <<- 0

	
  if (isTRUE(identify==FALSE))
    identify = NULL
    
  if (!is.null(identify)) {
    ivar = identify
       
    lab = rep("",m$nA)
    trellis.focus("panel", 1, 1)
    
    a.vec = NULL
    dev.new = NULL
    
    
    while(TRUE) {
      a = panel.identify(n=1,labels = lab)
      if (NROW(a)<1) break;
      ltext(m$action.val.mat[a,1],m$action.val.mat[a,2],labels=NROW(a.vec)+1,col=col.a[NROW(a.vec)+1])
      dev.act = dev.cur()
      if (is.null(dev.new) | NROW(dev.list()) <2) {
        win.graph()
        dev.new = dev.cur()
      } else {
        dev.set(dev.new)
      }
      setSimplexParmGLPK(IT_LIM,iter.lim)
      glp_std_basis(m$lp);
      
      if (ivar == "Ue") {
        ret = get.Ue.L.a(m,a, lp=m$lp, do.plot=FALSE)
      } else if (ivar == "v1") {
        ret = get.vi.L.a(m,a, i=1,lp=m$lp, do.plot=FALSE)
      } else if (ivar == "v2") {
        ret = get.vi.L.a(m,a, i=2,lp=m$lp, do.plot=FALSE)
      } else {
        stop("identify must be set to NULL, 'Ue', 'v1' or 'v2')")
      }        
  
         
      flush.console()
      mat = ret$LY.mat[,1:2,drop=FALSE]
      colnames(mat)=c("L",ivar)
      m$LUe[[a]] = mat

      xli[2] <<- max(c(mat[,"L"],xli[2]))
      yli[1] <<- min(c(mat[,ivar],yli[1]))
      yli[2] <<- max(c(mat[,ivar],yli[2]))
      
      a.vec = c(a.vec,a)
      
      plot(0,0,type="n",main=ivar,ylab=ivar,xlab="L",xlim=c(xli[1],xli[2]+diff(xli)*0.3), ylim=yli )
      
      for (i in 1:NROW(a.vec)) {
        lines(m$LUe[[a.vec[i]]][,1],m$LUe[[a.vec[i]]][,2],col=col.a[i],lty=1)
        points(m$LUe[[a.vec[i]]][,1],m$LUe[[a.vec[i]]][,2],col=col.a[i])
      }
      for (i in NROW(a.vec):1) {
        lines(m$LUe[[a.vec[i]]],col=col.a[i],lty=2)
      }
    
      legend("topright",legend=m$lab.a[a.vec],fill=col.a[1:NROW(a.vec)])
      
      dev.set(dev.act)
      bringToTop(which = dev.act, stay = FALSE)
    }
  }
  invisible(return(m))
}  


panel.levelplot.repgame = function(x, y,z,at,m,region=TRUE, contour=FALSE,
           col.nash="purple",col.max.G="white",col.br1 = "green", col.br2 = "greenyellow",
           col.ae = hsv(2/3,0.6,0.6), col.a1 = hsv(0,0.6,0.6), br.as.lines=FALSE, ...) {
  panel.levelplot(x,y,z,at=at,region=region,contour=contour,...)
  
 
  if (!is.null(m$opt.mat)) {
    
    #col = hsv(2/3,1/2,1/2,alpha=0.5)
    
    if (!is.na(col.ae)) {
      rows = m$opt.mat[,"opt"]==1
      ae.ind = m$opt.mat[rows,"ae"]

      lpoints(m$action.val.mat[ae.ind,1],m$action.val.mat[ae.ind,2],col=col.ae,pch=19)
      llines(m$action.val.mat[ae.ind,1],m$action.val.mat[ae.ind,2],col=col.ae,pch=19)
    }
    
    if (!is.na(col.a1)) {
      rows = m$opt.mat[,"opt"]==1
      a1.ind = m$opt.mat[rows,"a1"]
      lpoints(m$action.val.mat[a1.ind,1],m$action.val.mat[a1.ind,2],col=col.a1,pch=19)
      llines(m$action.val.mat[a1.ind,1],m$action.val.mat[a1.ind,2],col=col.a1,pch=19) 
    }   
  }
   
  # Best reply of player 1
  if (!is.na(col.br1)) {
    rows = m$c[,1] == m$g[,1]
    lpoints(m$action.val.mat[rows,1],m$action.val.mat[rows,2],col=col.br1,pch=1)
    #llines(m$action.val.mat[rows,1],m$action.val.mat[rows,2],col=col.br1,lwd=2,pch=3)
  }
    
  # Best reply of player 2
  if (!is.na(col.br2)) {
    rows = m$c[,2] == m$g[,2]
    lpoints(m$action.val.mat[rows,1],m$action.val.mat[rows,2],col=col.br2,pch=1)
  }
  
  
  # Stage game nash equilibrium  
  if (!is.na(col.nash)) {
    lpoints(m$action.val.mat[m$nash,1],m$action.val.mat[m$nash,2],col=col.nash,pch=19)
    lpoints(m$action.val.mat[m$nash,1],m$action.val.mat[m$nash,2],col=col.nash,pch=4)
    lpoints(m$action.val.mat[m$nash,1],m$action.val.mat[m$nash,2],col=col.nash,pch=3)
  }
}



# Draws convex hull of stage game payoffs and payoff set of the repeated game with side payments
plot.payoff.set = function(m, delta, main=NULL, plot.all.points = NROW(m$G)<=16,plot.a.lab = 1,u1.min=NULL,u1.max=NULL,u2.min=NULL,u2.max=NULL) {
	
	#restore.point("plot.payoff.set")
	
	if (m$n > 2) {
		warning("plot.payoff.set only works for two-player games")
		return()
	}
	if (is.null(main)) {
		main = paste(m$name," (Delta=", delta,")",sep="")
	}
	
	
	v1 = get.payoff.at.delta(m,delta,"v1")
	v2 = get.payoff.at.delta(m,delta,"v2")
	Ue = get.payoff.at.delta(m,delta,"Ue")

	u.minmax = c(min(m$c[,1]),min(m$c[,2]))	
	
	
	xlim=range(c(m$g[,1],Ue-v2)); 
	ylim=range(c(m$g[,2],Ue-v1));
	xlim[1] = u.minmax[1]; ylim[1] = u.minmax[2];
	ylim[1] = ylim[1]- 0.1 * diff(ylim); # Add 10% at the bottom to draw text labels
	xlim[1] = xlim[1]- 0.05 * diff(xlim); # Add 5% to the left to draw text labels
	xlim[2] = xlim[2]+ 0.04 * diff(xlim); # Add 4% to the right to draw text labels
	
	if (!is.null(u1.min)) xlim[1]=u1.min;
	if (!is.null(u2.min)) ylim[1]=u2.min;
	if (!is.null(u1.max)) xlim[2]=u1.max;
	if (!is.null(u2.max)) ylim[2]=u2.max;
	
	plot(m$g[,1],m$g[,2],type="n",pch=19,main=main,xlab="u1",ylab="u2",xlim=xlim,ylim=ylim)
	# Draw convex hull of stage game payoffs
	ch = chull(m$g)

	if (plot.all.points) {
		points.ind = 1:m$nA
	}	else {
		points.ind = ch
	}
	
	
	polygon(x=c(v1,v1,Ue-v2),y=c(v2,Ue-v1,v2), col="lightblue",border="blue",lwd=2)
	polygon(m$g[ch,1],m$g[ch,2], border="black")
	
		
	points(m$g[points.ind,1],m$g[points.ind,2],type="p",pch=19)
	
	points(m$g[m$nash,1],m$g[m$nash,2],type="p",pch=19,col="blue")

	points(u.minmax[1],u.minmax[2],type="p",pch=19,col="red")
	
	abline(h=u.minmax[1],lty=2)
	abline(v=u.minmax[2],lty=2)
	#polygon(x=c(v1,v1,Ue),y=c(v2,Ue,v2),border="blue",lwd=2)
	if (plot.a.lab) {
		if (!plot.all.points) {
			ind = which.labels.to.draw(x=m$g[ch,1],y=m$g[ch,2],xlim,ylim,min.xd = 0.07,min.yd=0.07)
			ind = ch[ind]
		} else {
			ind = points.ind
		}
		text(m$g[ind,],paste("(",m$lab.a[ind],")",sep=""),pos=1,cex=0.8)
	}

}


levelplot.Ue.L = function(m,L,main=paste(m$name, " Ue (L=",L,")", sep=""),cuts=100,
                          panel = panel.levelplot.repgame,...) {
  Ue.L = all.a.get.Ue.L(m,L)
  
  sk.levelplot(grid.xyz = data.frame(m$action.val.mat, Ue.L),
             main=main,cuts=cuts, m=m,
             panel = panel,...)
}


which.labels.to.draw = function(x,y,xlim,ylim,min.xd=0.1,min.yd=0.1) {
	
	#restore.point("which.labels.to.draw")

	if (length(xlim)==2) xlim=xlim[2]-xlim[1];
	if (length(ylim)==2) ylim=ylim[2]-ylim[1];
	
	ind = 1:NROW(x)
	while(NROW(x)>1) {
		ind0 = 1:(NROW(x)-1)
		ind1 = ind0+1
		nearer = c(FALSE, abs(x[ind0]-x[ind1])/xlim < min.xd & abs(y[ind0]-y[ind1]) < min.yd*ylim)
		del = rep(FALSE,NROW(x))
		del[-1] = nearer[-1] & (!nearer[-NROW(x)] | ((2:NROW(x) %% 2) == 0))
		if (sum(del)==0) break;
		x = x[!del]; y=y[!del]; ind = ind[!del];
	}
	return(ind)
}

		
		
	
	

draw.line.with.point.and.slope = function(X,Y,slope,...) {
	a = Y-X*slope
	abline(a=a,b=slope,...)
}

plot.LY.mat = function(LY.mat,LY.bound.mat=NULL,xlab="L",...) {
	
#	assign.list(VAR.STORE[["plot.LY.mat"]]$var)		
	
	plot(LY.mat[,"L"],LY.mat[,"Y"], ylim=range(LY.mat[,"Y"],LY.bound.mat[,"Y"],na.rm=TRUE),xlab=xlab,...)
	lines(LY.mat[,"L"],LY.mat[,"Y"], lty=2)
	kinks = LY.mat[,"kink"] == 1
	points(LY.mat[kinks,"L"],LY.mat[kinks,"Y"], pch=19, col="blue", lty=2)
	if (!is.null(LY.bound.mat)) {
		points(LY.bound.mat[,"L"],LY.bound.mat[,"Y"], pch=19, col="green", lty=2)	
	}
}	



#' Shows a plot that illustrates solution of a repeated game

#' @param m solved model
#' @param xvar variable to show on the x-Axis (default="delta", discount factor), should be contained in m$opt.mat or in mat1 and mat2
#' @param yvar variable to show on the y-Axis (default="Ue", joint payoffs), should be contained in m$opt.mat or in mat1 and mat2
#' @param identify if true can click on plot to get additional information on corresponding payoffs and optimal action plans corresponding to that x-Value (in particular for xvar="delta")
#' @param legend.pos position of legend, set NULL for no legend
#' @param ... further parameters for plot
plot.repgame = function(m,xvar="delta",yvar=NULL,main=m$name,
                        delta.seq = NULL, r.seq=NULL, 
                        add.crit.delta=TRUE,add.last.L = TRUE,identify=NULL,...) {
	restore.point("plot.repgame")
	#restore.point("plot.repgame")

	                        	                        
	if (is.null(m$opt.mat)) {
		warning("plot.repgame: Model not yet solved. Cannot plot.")
		return();
	}

	if (xvar=="delta") {
		#if (is.null(delta.seq) & (m$sol.type == "imp" | !is.null(identify))) {
  	if (is.null(delta.seq)) {
			delta.seq = seq(0,1,by=0.001)
		}
		mat = get.mat.of.delta(m=m,delta=delta.seq,add.crit.delta=add.crit.delta)
	}	else if (xvar=="r") {
		if (is.null(r.seq)) {
			delta.seq = seq(0,1,by=0.001)
		} else {
			delta.seq = 1 / 1+r.seq
		}
		mat = get.mat.of.delta(m=delta,delta=delta.seq,add.crit.delta=add.crit.delta)
	} else {
		mat = m$opt.mat
		if (add.last.L) {
			mat = rbind(mat,mat[NROW(mat),])
			mat[NROW(mat),"L"] = mat[NROW(mat),"L"]* 1.15
		  mat[NROW(mat),"r"] = (mat[NROW(mat),"Ue"]-mat[NROW(mat),"V"]) / mat[NROW(mat),"L"]
			mat[NROW(mat),"delta"] = 1 / (1+mat[NROW(mat),"r"])
			mat[NROW(mat),"opt"] = 0
		}
	}
	plot.mat.repgame(m,mat=mat,xvar=xvar,yvar=yvar,main=main,identify=identify,...)
}

plot.mat.repgame = function(m=NULL,mat, legend.pos ="topleft",xvar = "delta", yvar=NULL,
														xlab = NULL, ylab=NULL,xlim=NULL, col=NULL,
                            ylim = NULL,lwd=1, add = TRUE, main=m$name,
                            draw.U.nash=TRUE, draw.U.max = ("Ue" %in% yvar), identify = NULL) {
	
	#restore.point("plot.mat.repgame")

	# Initialize default values for ranges, labels and colors

	if (is.null(yvar))
  	yvar = c("Ue","V")
	      	                            
	if (is.null(xlim))
		xlim = range(mat[,xvar])

	if (m$sol.type == "pmx") {
	  draw.U.max = FALSE
    draw.U.nash = FALSE
  }		
  
	if (is.null(ylim)) {
		ylim = range(mat[,yvar])
		if (m$sol.type != "pmx") {
  		if ("Ue" %in% yvar) 
  			ylim[2] = max(m$G)
  		if ("V" %in% yvar)
  			ylim[1] = m$V.min.max
  		if (NROW(yvar)==1 & yvar[1]=="Ue" & NROW(m$nash)>=1)
    		ylim= c(max(m$G[m$nash]),max(m$G))
  	}
  	if (!is.null(legend.pos)) {
      if (legend.pos=="topleft" | legend.pos == "topright") {
        if (draw.U.max) {
          ylim[2] = ylim[2]+0.25*diff(ylim)
        } else {
          ylim[2] = ylim[2]+0.2*diff(ylim)
        }
      }
    }               		
	}
				
	if (is.null(xlab)) {
		xlab = REPGAME.XLAB.DEFAULT[xvar]
		if (is.na(xlab)) {
		  xlab = xvar
	  }
  }
		
	if (is.null(ylab)) {
		if (!is.null(legend.pos)) {
			ylab=""
		} else {
			ylab = paste(yvar,collapse=" ")
		}
	}
	
	if (is.null(col)) {
		col = rep("black",NROW(yvar))
		col[match("Ue",yvar,nomatch=NULL)] = "blue"
		col[match("V",yvar,nomatch=NULL)] = "red"
	}

	if (add) {
		plot(0,0, type="n",xlab=xlab,ylab=ylab, main=main,
	  	   xlim=xlim,ylim=ylim, col="blue")
 	}
 	
 	#draw Nash and max payoff
 	if (draw.U.nash)
   	if (NROW(m$nash>1)) 
     	abline(h=max(m$G[m$nash]),lty=2,col="grey")
 	if (draw.U.max)
   	abline(h=max(m$G),lty=2,col="grey")
     	
 	
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
		legend(legend.pos, legend=yvar, fill=col)
	}
	
  if (isTRUE(identify==FALSE))
    identify = NULL
	if (!is.null(identify)) {
  	x = mat[,xvar]
  	while(TRUE) {
    	# Identify x value
  	  ind = identify(x, y = rep(mean(ylim),NROW(x)),pos = FALSE,
                     n = 1, plot = FALSE, atpen = FALSE, offset = 0.5,
                     tolerance = 50,labels=rownames(mat))
      if (NROW(ind)<1)
        break()
      #rownames(both)=c(m1.name,m2.name)                    
    	print(paste(xvar, "=", x[ind]))
    	print(mat[ind,,drop=FALSE])
    	
    	flush.console()

    	abline(v = x[ind], col="grey", lty=2)
  	}
	}	
}


REPGAME.XLAB.DEFAULT = c("Discount Factor", "Discount Rate","Liquidity L")
names(REPGAME.XLAB.DEFAULT)=c("delta","r","L")


#' Shows a plot that compares several models in one plot

#' @param m.li a list containing several models
#' @param xvar variable to show on the x-Axis (default="delta", discount factor), should be contained in m$opt.mat or in mat1 and mat2
#' @param yvar variable to show on the y-Axis (default="Ue", joint payoffs), should be contained in m$opt.mat or in mat1 and mat2
#' @param x.seq values on the x-Axis
#' @param x.min minimum value on x-Axis
#' @param x.max maximum value on x-Axis
#' @param identify if true can click on plot to get additional information on corresponding payoffs and optimal action plans corresponding to that x-Value (in particular for xvar="delta")
#' @param legend.pos position of legend, set NULL for no legend

plot.models = function(m.li,mat.li=NULL,names = unlist(sublist(m.li,"name")), xvar = "delta",yvar=c("Ue"),
                       x.seq=NULL, x.min=NULL,x.max = NULL,
                       col=NULL,legend.pos ="topright",lwd=1, eps =  REP.CNT$TOL,
                       identify=NULL,main=NULL,ylab=NULL,xlab=NULL) {
                             
  restore.point("plot.compare.models")
	#restore.point("plot.compare.models")                             
  
	
  if (is.null(mat.li)) {
    mat.li=list()           
    if (xvar == "delta") {
    	if (is.null(x.seq))
  	  	x.seq = seq(0,1,by=0.00001)
      
  	  if (!is.null(x.min)) x.seq = x.seq[x.seq>=x.min]
  	  if (!is.null(x.max)) x.seq = x.seq[x.seq<=x.max]
  	  
  	  for (i in 1:length(m.li)) {
  	    mat.li[[i]] = get.mat.of.delta(m=m.li[[i]],delta=x.seq,add.crit.delta=FALSE, pm.keep.sparse = FALSE)
	    }
    } else if (xvar == "L") {
      opt.mat.li = sublist(m.li,opt.mat) 
      joined.opt.mat = rbind.list(opt.mat.li, only.common.cols = TRUE)

    	if (is.null(x.seq)) {	
      	L.min = 0; L.max = max(joined.opt.mat[,"L"])
      	L = sort(unique(c(seq(0,L.max,length=100),
      	              joined.mat[,"L"],joined.mat[joined.mat[,"opt"]==0,"L"]-eps)))
        x.seq = L        
        if (!is.null(x.min)) x.seq = x.seq[x.seq>=x.min]
    	  if (!is.null(x.max)) x.seq = x.seq[x.seq<=x.max]
  	  }
      for (i in 1:length(m.li))
        mat.li[[i]] = get.mat.of.L(m.li=m.li[[1]],L=x.seq,add.crit.L=FALSE)
      
    } 
  } else {
    x.seq = mat1[,"xvar"]
  }
  joined.mat = rbind.list(mat.li,only.common.cols = TRUE)
  
  ylim = range(joined.mat[,yvar])
  if (diff(ylim)==0) {
    ylim[2] = ylim[1]+1
  } else {
    if (!is.null(legend.pos)) {
      if (legend.pos == "topleft" | legend.pos == "topright") {
        ylim[2] = ylim[2]+0.25*diff(ylim)
      }
    }
  }
    
  xlim = range(x.seq)
  if (is.null(ylab))
    ylab = paste(yvar,collapse=" ")
  if (is.null(main)) 
    main=paste(names, sep=" vs ")
  if (is.null(xlab))
    xlab=REPGAME.XLAB.DEFAULT[xvar]
    
  if (is.null(col)) {
    col = c("black","blue","red","green","orange","purple","grey","yellow")[length(m.li)]
	}

	draw.plot = function () {
   	plot(0,0, type="n",xlab=xlab,ylab=ylab, main=main,
  	  	   xlim=xlim,ylim=ylim, col="black")

  	for (i in 1:length(m.li)) {
      lines(mat.li[[i]][,xvar],mat.li[[i]][,yvar],col=col[i],lty=1)
    }
    for (i in (length(m.li)-1):1) {
      lines(mat.li[[i]][,xvar],mat.li[[i]][,yvar],col=col[i],lty=2)
    }
    
    if (!is.null(legend.pos)) {
  		legend(legend.pos, legend=paste(yvar,names),
  		       fill=col)
		}
	}
	draw.plot()
}  

#' Shows a plot that compares two models

#' @param m1 first model
#' @param m2 second model
#' @param xvar variable to show on the x-Axis (default="delta", discount factor), should be contained in m$opt.mat or in mat1 and mat2
#' @param yvar variable to show on the y-Axis (default="Ue", joint payoffs), should be contained in m$opt.mat or in mat1 and mat2
#' @param x.seq values on the x-Axis
#' @param x.min minimum value on x-Axis
#' @param x.max maximum value on x-Axis
#' @param identify if true can click on plot to get additional information on corresponding payoffs and optimal action plans corresponding to that x-Value (in particular for xvar="delta")
#' @param legend.pos position of legend, set NULL for no legend
plot.compare.models = function(m1,m2,xvar = "delta",yvar=c("Ue"),x.seq=NULL, x.min=NULL,x.max = NULL,
                           col=NULL,legend.pos ="topright",lwd=1, eps =  REP.CNT$TOL,
                           m1.name =  "m1", m2.name="m2", mat1=NULL,mat2 = NULL,
                           show.diff=FALSE, identify=NULL,main=NULL,ylab=NULL,xlab=NULL) {
                             
  restore.point("plot.compare.models")
	#restore.point("plot.compare.models")                             
  
  if (is.null(mat1)) {               
    if (xvar == "delta") {
    	if (is.null(x.seq))
  	  	x.seq = seq(0,1,by=0.00001)
      
  	  if (!is.null(x.min)) x.seq = x.seq[x.seq>=x.min]
  	  if (!is.null(x.max)) x.seq = x.seq[x.seq<=x.max]
  	  
    	mat1 = get.mat.of.delta(m=m1,delta=x.seq,add.crit.delta=FALSE, pm.keep.sparse = FALSE)
  	  mat2 = get.mat.of.delta(m=m2,delta=x.seq,add.crit.delta=FALSE, pm.keep.sparse = FALSE)
    } else if (xvar == "L") {
    	if (is.null(x.seq)) {
      	L.min = 0; L.max = max(m1$opt.mat[,"L"],m2$opt.mat[,"L"])
      	L = sort(unique(c(seq(0,L.max,length=100),
      	                  m1$opt.mat[,"L"],m1$opt.mat[m1$opt.mat[,"opt"]==0,"L"]-eps,
      	                  m2$opt.mat[,"L"],m2$opt.mat[m2$opt.mat[,"opt"]==0,"L"]-eps)))
  
        x.seq = L
        
        if (!is.null(x.min)) x.seq = x.seq[x.seq>=x.min]
    	  if (!is.null(x.max)) x.seq = x.seq[x.seq<=x.max]
  
      	mat1 = get.mat.of.L(m=m1,L=x.seq,add.crit.L=FALSE)
  	    mat2 = get.mat.of.L(m=m2,L=x.seq,add.crit.L=FALSE)
      }
    } 
  } else {
    x.seq = mat1[,"xvar"]
  }


	if (show.diff) {
  	mat = mat1    
  	mat[,yvar] = mat1[,yvar]-mat2[,yvar]
    ylim = range(mat[,yvar])
  } else {
    ylim = range(c(mat1[,yvar],mat2[,yvar]))
  }
  
  
  if (diff(ylim)==0) {
    ylim[2] = ylim[1]+1
  } else {
    if (!is.null(legend.pos)) {
      if (legend.pos == "topleft" | legend.pos == "topright") {
        ylim[2] = ylim[2]+0.25*diff(ylim)
      }
    }
  }
    
  xlim = range(x.seq)
  if (is.null(ylab))
    ylab = paste(yvar,collapse=" ")
  if (is.null(main)) 
    main=paste(m1.name,"vs",m2.name)
  if (is.null(xlab))
    xlab=REPGAME.XLAB.DEFAULT[xvar]
    
  if (is.null(col)) {
		col = rep("black",NROW(yvar))
		col[match("Ue",yvar,nomatch=NULL)] = "blue"
		col[match("V",yvar,nomatch=NULL)] = "red"
	}

	draw.plot = function () {
   	plot(0,0, type="n",xlab=xlab,ylab=ylab, main=main,
  	  	   xlim=xlim,ylim=ylim, col="black")
  	  	   
  	if (show.diff) {
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
    		legend(legend.pos, legend=yvar, fill=col)
  	  }

    } else {
      lines(mat1[,xvar],mat1[,yvar],col="blue",lty=1)
      lines(mat2[,xvar],mat2[,yvar],col="red",lty=1)
      lines(mat1[,xvar],mat1[,yvar],col="blue",lty=2)
      
      if (!is.null(legend.pos)) {
    		legend(legend.pos, legend=paste(yvar,c(m1.name,m2.name)),
    		       fill=c("blue","red"))
  	  }
    }
  
	}
	draw.plot()

	draw.two.Ue.L.a = function(m1,m2,a.vec,L,Ue1,Ue2) {
  	restore.point("draw.two.Ue.L.a")
	  #restore.point("draw.two.Ue.L.a")
  	
	  dev.act = dev.cur()
    win.graph()
   
    
    ret = list()
    ret[[1]] = get.Ue.L.a(m1,a.vec[1], lp=m1$lp, do.plot=FALSE,Y.min=-Inf)
    flush.console()
    ret[[2]] = get.Ue.L.a(m2,a.vec[2], lp=m2$lp, do.plot=FALSE,Y.min=-Inf)
    flush.console()
    
    mat1 = ret[[1]]$LY.mat[,1:2,drop=FALSE]; 
    mat2 = ret[[2]]$LY.mat[,1:2,drop=FALSE]; 
    colnames(mat1)=colnames(mat2) = c("L","Ue")
    
    xli = c(0,max(c(mat1[,"L"],mat2[,"L"])))
    yli = c(min(c(mat1[,"Ue"],mat2[,"Ue"])),max(c(mat1[,"Ue"],mat2[,"Ue"])))
    
    main = paste("Ue a(m1)=",m1$lab.a[a.vec[1]]," a(m2)=",m2$lab.a[a.vec[2]])
    plot(0,0,type="n",main=main,ylab="Ue",xlab="L",xlim=c(xli[1],xli[2]+diff(xli)*0.3), ylim=yli )

    abline(v=L,col="grey",lty=2)
    abline(h=Ue1,col="lightblue",lty=2)
    abline(h=Ue2,col="orange",lty=2)
    lines(mat1[,"L"],mat1[,"Ue"],col="blue",lty=1)
    lines(mat2[,"L"],mat2[,"Ue"],col="red",lty=1)
    lines(mat1[,"L"],mat1[,"Ue"],col="blue",lty=2)
      
    dev.set(dev.act)
    bringToTop(which = dev.act, stay = FALSE)
  }	
	
  
	
  if (isTRUE(identify==FALSE))
    identify = NULL
    
  if (!is.null(identify)) {
    x = x.seq
  	while(TRUE) {
  	  ind = identify(x, y = rep(mean(ylim),NROW(x)),pos = FALSE,
                     n = 1, plot = FALSE, atpen = FALSE, offset = 0.5,
                     tolerance = 50)
      if (NROW(ind)<1)
        break()
       
      if (m1$sol.type == m2$sol.type & NCOL(mat1) == NCOL(mat2)) {
        both = rbind(mat1[ind,,drop=FALSE],mat2[ind,,drop=FALSE])
        rownames(both)=c(paste(m1.name,rownames(mat1)[ind]),
                         paste(m2.name,rownames(mat2)[ind]))   
        print(paste(xvar, "=", x[ind]))
      	print(both)
    	} else {
        print(paste(xvar, "=", x[ind]))
        print(m1.name)
        print(mat1[ind,,drop=FALSE])
        print(m2.name)
        print(mat2[ind,,drop=FALSE])
      }
        
    	flush.console()
    	if (m1$sol.type == "imp" & m2$sol.type=="imp") {
        draw.plot()
        draw.two.Ue.L.a(m1,m2,c(mat1[ind,"ae"],mat2[ind,"ae"]),L=x[ind],Ue1=mat1[ind,"Ue"],Ue2=mat2[ind,"Ue"])
  	  }
      abline(v = x[ind], col="grey", lty=2)
    	
    	flush.console()
  
  	}
	}
}
 	   
	

plot.Ue.opt = function(Ue.opt,xlim=NULL,ylim=NULL,xlab=colnames(Ue.opt)[1],ylab=colnames(Ue.opt)[2],...) {
	
	#assign.list(VAR.STORE[["plot.Ue.opt"]]$var)		
	
	X = Ue.opt[,1]
	Y = Ue.opt[,2]

	if (is.null(xlim)) {
		xlim = range(X,na.rm=TRUE)
	}
	if (is.null(xlim)) {
		ylim = range(Y,na.rm=TRUE)
	}			
	plot(X,Y,xlim=xlim,ylim=ylim, xlab=xlab,ylab=ylab,...)
	lines(Ue.opt,lty=2)
}
			
