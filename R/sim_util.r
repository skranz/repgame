# Utility functions to run simulations of some numerical model
# for different parameter combinations

# First generate a parameter table

para.block.to.matrix = function(block,para.def=NULL) { 
  para.def = copy.into.list(dest=para.def,source=block$def)
  para.def = unlist(para.def)
  val.dim = sapply(block$val,length)
  nrows = prod(val.dim)
  mat = matrix(para.def,nrows,length(para.def),byrow=TRUE)
  colnames(mat)=names(para.def)
  mat[,names(block$val)] = make.grid.matrix(x=block$val)
  return(mat)
}
                        
string.const.to.int.const = function(str,env=globalenv(),prefix=NULL) {
  if (is.list(str)) {
    for (i in 1:length(str)) {
      string.const.to.int.const(str[[i]],env=env,prefix=prefix)
    }
    return();
  }
  for (i in 1:length(str)) {
    na = str[i]
    na = paste(prefix,toupper(na),sep="")
    assign(na,i, envir=env)
  }
}

get.changed.cols = function(mat,row=NULL) {
  changed = NULL
  if (row >1) {
    changed = which(mat[row-1,] != mat[row,])
  }
  if (row<NROW(mat)) {
    changed = union(changed,which(mat[row-1,] != mat[row,]))
  }
  return(changed)
}

#COLL = 0
#string.const.to.int.const(list(c("COLL","MON","FB")))
#COLL


                        
new.sim = function(para.def,blocks, para.names = names(para.def),para.factors=NULL,file.name=NULL) {
  sim = ListEnv()
  copy.into.env(dest=sim,exclude="sim") # Copy all parameters into sim

  para.mat = NULL
  block.id = NULL
  for (i in 1:length(blocks)) {
    block.mat = para.block.to.matrix(blocks[[i]],para.def)
    para.mat = rbind(para.mat,block.mat)
    block.id = c(block.id,rep(i,NROW(block.mat)))
  }
  sim$para.mat = para.mat
  sim$run.mat = matrix(0,NROW(para.mat),4)
  colnames(sim$run.mat) = c("run.id","block","model","solved")
  sim$run.mat[,"block"] = block.id
  sim$run.mat[,"model"] = 1:NROW(sim$run.mat)
  
  sim$file.name = file.name
  class(sim)=c("sim","ListEnv")
  return(sim)
}


run.sim = function(sim, file.name = sim$result.file, append=TRUE,
                   load.data = TRUE, rows=NULL) {
  restore.point("run.sim")
  
  stopifnot(!is.null(file.name))
  
  sim$do.update = 0
  sim$run.id = as.numeric(Sys.time())
  sim$run.mat[,"run.id"] = sim$run.id
  sim$file.name = file.name
  
  sim$result.file = paste("res_",sim$file.name,sep="")
  sim$model.file =  paste("mod_",sim$file.name,sep="")
  
  num.s = NROW(sim$para.mat)
  
	prev.para = NULL
	sol.added = FALSE
	
	
  sim.par.s = function(s) { 
		restore.point("sim.par.s")
	
    para.vec = sim$para.mat[s,]
    para = as.list(sim$para.mat[s,])
    print("Now...")
    print(unlist(para))
    if (s==1) {
      sim$init(sim,para)
    } else {
      sim$update(sim,para,prev.para)
    }
		prev.para <<- para
    block.par = names(sim$blocks[[sim$run.mat[s,"block"]]]$val)
    sim$m$name = paste("(",s,"/",num.s,")",
                       paste(block.par, para.vec[block.par],sep="=",collapse=" "))
    
		
		sim$solve(sim,para)
    res = sim$result(sim,para)
    
    # No solution found
    if (is.null(res)) {
      col.names = (s==1) & ((!file.exists(sim$result.file)) | (!append))    
      runpara  = c(sim$run.mat[s,],sim$para.mat[s,])
      sim$run.mat[s,"solved"]=-1
      # Write entry into model file
      write.table(t(runpara), file = sim$model.file, append = append | s>1,
        quote = FALSE,sep = "\t",eol = "\n", na = "NA", dec = ".",
        row.names=FALSE, col.names = col.names)
      print(paste("No solution found for model",s))
      print(runpara)
      return()
    }

       
    sim$run.mat[s,"solved"]=1
    
    runpara  = c(sim$run.mat[s,],sim$para.mat[s,])
    pm = matrix(runpara,NROW(res),length(runpara),byrow=TRUE)
    colnames(pm)=c(colnames(sim$run.mat),colnames(sim$para.mat))
    mat = cbind(pm,res)
    # Report to output
    col.names = (!sol.added) & ((!file.exists(sim$result.file)) | (!append))
    
    runpara  = c(sim$run.mat[s,],sim$para.mat[s,])

    # Write results into a file
    if (!is.null(sim$result.file)) {
      if (!sol.added) {		
        write.table(mat, file = sim$result.file, append = append, quote = FALSE, 
            sep = "\t",eol = "\n", na = "NA", dec = ".", row.names=FALSE, col.names = col.names)
      } else {
        write.table(mat, file = sim$result.file, append = TRUE, quote = FALSE, 
            sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
      }

      # Write entry into model file
			col.names = (!sol.added) & ((!file.exists(sim$result.file)) | (!append))
      write.table(t(runpara), file = sim$model.file, append = append | s>1,
        quote = FALSE,sep = "\t",eol = "\n", na = "NA", dec = ".",
        row.names=FALSE, col.names = col.names)
				
			sol.added <<- TRUE
    }
    print("writen to file...")
  }

    
  if (is.null(rows)) {
    rows = 1:NROW(sim$para.mat)
  }
  for (s in rows) {
    sim.par.s(s)
  }
  if (load.data) {
    load.sim(sim=sim)
  }
}

# Makes only sense if solution has been modified, e.g. duplicated rows removed
save.sim = function(sim,file.name=sim$file.name, append = FALSE) {
  sim$file.name = file.name
  
  sim$result.file = paste("res_",sim$file.name,sep="")
  sim$model.file =  paste("mod_",sim$file.name,sep="")
  col.names = ((!file.exists(sim$result.file)) | (!append))
  write.table(sim$sol, file = sim$result.file, append = append, quote = FALSE, 
    sep = "\t",eol = "\n", na = "NA", dec = ".", row.names=FALSE, col.names = col.names)
  write.table(sim$mod.mat, file = sim$result.file, append = append, quote = FALSE, 
    sep = "\t",eol = "\n", na = "NA", dec = ".", row.names=FALSE, col.names = col.names)
  print("sim saved....")
  
}


load.sim = function(sim=NULL, file.name=sim$file.name, return.sim=is.null(sim),as.matrix=!TRUE) {
  
  #restore.point("load.sim")
  
  if (is.null(sim)) {
    return.sim =TRUE
    sim = ListEnv()
  }
  
  sim$file.name=file.name
  sim$result.file = paste("res_",sim$file.name,sep="")
  sim$model.file =  paste("mod_",sim$file.name,sep="")

  sol = read.table(sim$result.file, header = TRUE, sep = "\t", quote = "\"'",
           dec = ".")
  mod.mat = read.table(sim$model.file, header = TRUE, sep = "\t", quote = "\"'",
           dec = ".")

  if (as.matrix) {
    sol = as.matrix(sol)
    mod.mat = as.matrix(mod.mat)
  }
 
  sim$sol = sol
  sim$mod.mat = mod.mat
    
  sim$para.mat = sim$mod.mat[,-(1:4)]
  sim$run.mat = sim$mod.mat[,(1:4)]
  if (return.sim) 
    return(sim)
}

new.dyngame.sim = function(para.def, para.names = names(para.def), para.factors=NULL,
                           blocks, make.my.game, change.tau=TRUE, change.g=TRUE, no.g.para=NULL, no.tau.para=NULL) {
  
  sim = new.sim(para.def,blocks=list(block))                           
  sim$init = function(sim=NULL,para) {
    print("update game...")
    my.game = make.my.game(para=para)
    sim$make.my.game = make.my.game
    sim$m = init.game(my.game=my.game)
  }
  sim$update = function(sim,para,prev.para = NULL,update.type=NULL) {
    restore.point("sim.update")
    #restore.point("sim.update")
    print("update model...")  
    m = sim$m
    my.game = make.my.game(para=para)
		
		# Look which parameters have been changed and check whether g or tau have to be changed
		if (!is.null(prev.para) & ((!is.null(no.g.para)) | !(!is.null(no.tau.para)))) {
			change.para = names(para)[which(unlist(prev.para)!=unlist(para))]
			change.g = length(intersect(change.para, setdiff(names(para),no.g.para)))>0
			change.tau = length(intersect(change.para, setdiff(names(para),no.tau.para)))>0
		}
		
    set.g.and.tau(m,g.fun=my.game$g.fun,tau.fun=my.game$tau.fun,recalc.g=change.g,recalc.tau = change.tau)
    m$delta = para$delta
    m$integrated = my.game$integrated
    print(paste("m$delta:", m$delta))
  }  
  sim$solve = function(sim,para) {
    sim$ms = solve.game(sim$m)
  }
  sim$result = function(sim,para) {
		restore.point("sim$result")
		#restore.point("sim$result")
		
		ms= sim$ms
    if (!ms$sol.exists) {
      return(NULL)
    }
    mat = cbind(ms$sol.mat,ms$extra.sol.cur,ms$extra.sol.ad,sim$m$xv.mat)
    return(mat)
  }
  sim$report = function(sim,para) {
    print("")
    print(unlist(para))
    print("solved")
    print("=========================================================================")
    print("=========================================================================")  
  }
  return(sim)
}

# Removes those solutions that have been duplicated
remove.duplicated.sim = function(sim,save=FALSE) { 
  dup.mod = duplicated(sim$mod.mat[,-(1:3)])
  dup.sol = duplicated(sim$sol[,-(1:3)])
  has.dup = any(c(dup.mod,dup.sol))
  
  if (has.dup) {
    if (any(sim$mod.mat[,"solved"]!=1)) {
      ord = order(sim$sol[,"solved"],decreasing=TRUE)
      sim$sol = sim$sol[ord,]
      dup.sol = dup.sol[ord]
      ord = order(sim$mod.mat[,"solved"],decreasing=TRUE)
      sim$mod.mat = sim$mod.mat[ord,]
      dup.mod = dup.mod[ord]
    }

    sim$mod.mat = sim$mod.mat[!dup.mod,]
    print(paste(sum(dup.mod), " duplicated models removed"))
    
    sim$sol = sim$sol[!dup.sol,]
    print(paste(sum(dup.sol), " duplicated solution rows removed"))

    if (save) {
      save.sim(sim)
    }  
    sim$para.mat = sim$mod.mat[,-(1:4)]
    sim$run.mat = sim$mod.mat[,(1:4)]
  }
}
 
sol.levelplot = function(sol,para,xyz,name="", arrows = TRUE,
                           xlab="x1", ylab="x2",main="",xrange=NULL,yrange=NULL,
                           clip.x=NULL,clip.y=NULL,...) {
  
  #restore.point("sol.levelplot")

  if (!is.list(para)) para = as.list(para)
  var =xyz
  para.pick = setdiff(names(para),var)
  rows = matrix.get.rows(sol,para[para.pick])
  mat = sol[rows,]
  
  if (!is.null(xrange)) {
    x.shown = which(mat[,var[1]]>=xrange[1] & mat[,var[1]]<=xrange[2])
    mat = mat[x.shown,]
  }
  if (!is.null(yrange)) {
    y.shown = which(mat[,var[2]]>=yrange[1] & mat[,var[2]]<=yrange[2])
    mat = mat[y.shown,]
  }  
 
  sk.levelplot(grid.xyz = mat[,var],xlab=xlab,ylab=ylab,main=main,...)
}


	
draw.sim.GUI = function(df, para = colnames(df), para.unique = NULL,max.option=4) {
	restore.point("draw.sim.GUI")
	#restore.point("draw.sim.GUI")
	
	
	var = colnames(df)
	require(gWidgets)
	require(gWidgetstcltk)
	options(guiToolkit="tcltk")

	win  <<-  gwindow("Sim Explorer",  visible=FALSE, width=150, height=80)
	nb <<- gnotebook(container=win)
			
	# Groups for the different pages
	#group.para <<- glayout(container=nb,horizontal=FALSE)
	group.main <<- ggroup(container=nb,horizontal=!FALSE)
	
	group.para <<- ggroup(container=group.main,horizontal=FALSE)
	group.plot <<- ggroup(container=group.main,horizontal=FALSE)

	if (is.null(para.unique)) {
		para.unique = list()
		for (i in 1:NROW(para)) {
			para.unique[[para[i]]] = unique(df[,para[i] ])
		}
	}
	
	# Sort such that parameters without variation are at the bottom
	para.length = sapply(para.unique,length)
	ord.para = order((para.length==1))
	para = para[ord.para]
	para.unique = para.unique[ord.para]
	para.length = para.length[ord.para]
	
	np = NROW(para)
	use.slider = sapply(1:np,
		function(i) {
			val = para.unique[[i]]
			if (NROW(val)>1) {
				if (all(diff(val)==diff(val[1:2]))) {
					return(diff(val[1:2]))
				}
			}
			return(0)
		}
	)
	

	handler = function(h,...) {
		ret = read.gui()
		do.plot(ret)
	}
	handler.para  = handler
	
	

	
	input.para.names = para
	
	label.para = list()
	input.para = list()
	cont = group.para
	handler = handler.para

	for (i in 1:NROW(para)) {

		val  = para.unique[[i]]		
		
		if (NROW(val)==1) {
			text = paste(para[i],": ",val[1],sep="")
		} else {
			text = paste(para[i],": ",sep="")
		}
		
		
		label.para[[para[i]]] = glabel(text = text, handler = NULL,container=cont)
		#cont[i,1,expand=!TRUE] = label.para[[i]]
		if (NROW(val)>max.option & use.slider[i]!=0) {
			input.para[[para[i]]] = gspinbutton(from=val[1],to=val[NROW(val)],by=val[2]-val[1], container=cont,handler=handler.para)
		} else if (NROW(val)>max.option) {
			input.para[[para[i]]] = gdroplist(items=val, selected = 1, editable = FALSE, container=cont,handler=handler.para)
		} else if (NROW(val)>=2) {
			input.para[[para[i]]] = gradio(items=val, selected = 1, horizontal = TRUE, container=cont,handler=handler.para)
		} else {
			#input.para[[i]] = gradio(items=val, selected = 1, horizontal = TRUE, container=cont,handler=handler.para)
		}
		#cont[i,2,expand=!TRUE] = input.para[[i]]	
	}
 
 
	# Plot window
	
	plot.types = c("image","plot")
	
	input.plot.names = c("plot.type","x1","x2","val")
	
	input.plot = list();
	label.plot = list();
	
	
	cont = group.plot;
	
	name = "plot.type";
	label.plot[[name]] = glabel(text = "plot.type", handler = NULL,container=cont)
	input.plot[[name]] = gdroplist(items=plot.types, selected = 1, editable = TRUE, container=	cont,handler=handler)

	name = "x1";	
	sel = which(para=="xv1")
	if (length(sel)<1) sel = 1
	label.plot[[name]] = glabel(text = "x (x1)", handler = NULL,container=cont)
	input.plot[[name]] = gdroplist(items=c("",para), selected = sel+1, editable = TRUE, container=	cont,handler=handler)

	
	name = "x2";	
	sel = which(para=="xv2")
	if (length(sel)<1) sel = 1
	label.plot[[name]] = glabel(text = "series (x2)", handler = NULL,container=cont)
	input.plot[[name]] = gdroplist(items=c("",para), selected = sel+1, editable = TRUE, container=	cont,handler=handler)

	name = "val";	
	label.plot[[name]] = glabel(text = "value", handler = NULL,container=cont)
	sel = which(var=="U")
	if (length(sel)<1) sel = 1
	input.plot[[name]] = gdroplist(items=var, selected = sel, editable = TRUE, container=	cont,handler=handler)

#	input.screen.x = gdroplist(items=para, selected = 1, editable = TRUE, container=		cont,handler=handler)
 
	#input.screen.y = gdroplist(items=c("",para), selected = 1, editable = TRUE, container=	cont,handler=handler)

	restore.point("draw.sim.GUI.later")
	#restore.point("draw.sim.GUI.later")
	
	read.gui = function() {
		#restore.point("draw.sim.GUI.later");
		
		# Values of inputs
		opt = list()
		
		opt$para = lapply(1:NROW(para), function(i) 
		{
			if (para.length[i]>1) {
				return(svalue(input.para[[para[i]]]))
			} else {
				return(para.unique[[para[i]]][1])
			}
		}
		)
		names(opt$para) = para
		opt$plot = lapply(names(input.plot), function(p) svalue(input.plot[[p]]))
		names(opt$plot) = names(input.plot)	
		return(opt)
	}

	do.plot = function(opt,sol=df) {
		restore.point("do.plot")
		#restore.point("draw.sim.GUI.later");restore.point("do.plot");
		
		type = opt$plot$plot.type
		if (type!="image") {
			warning("Only image plots implemented so far")
			return()
		}
		x1 = opt$plot$x1
		x2 = opt$plot$x2
		if (!(x1 %in% para)) return()
		if (!(x2 %in% para)) return()
		val = opt$plot$val
		if (!(val %in% var)) return()
	
	  sol.levelplot(sol=sol,para=opt$para,c(x1,x2,val),xlab=x1,ylab=x2,main=val,arrows=FALSE)
	}
	
	visible(win) <- TRUE
}

#para=list(x.max=30, x.learn=15, kappa=20, rho=0.5,gamma=0.1,sigma=0.1,delta=0.5,p.min=0,p.max=25,p.step=1,tau.sparse = TRUE, pred.pricing = !TRUE)
#para = c(names(para),"xv1","xv2")
#draw.sim.GUI(sim$sol,para=para)





