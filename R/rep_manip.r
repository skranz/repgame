
	make.game.with.signal.manipulation = function(m,pi.man,cost) {
  	restore.point("make.game.with.signal.manipulation")
  	# restore.point("make.game.with.signal.manipulation")
 
  	
		if (m$n>2) {stop("make.manipulation.game only implemented for 2 players")}
	
		# Number of actions in orignal game
		ona = c(NROW(m$g1),NCOL(m$g1))
		og1 = m$g1; og2 = m$g2
		ny = m$ny
		na = ona*(ny+1)


		# Mappings between new and old game actions		
		lab.a1 = make.grid.matrix(x=list(c("",m$lab.y),m$lab.ai[[1]]))
		lab.a1 = paste(lab.a1[,2],".",lab.a1[,1],sep="")
		lab.a2 = make.grid.matrix(x=list(c("",m$lab.y),m$lab.ai[[2]]))
		lab.a2 = paste(lab.a2[,2],".",lab.a2[,1],sep="")		
		
		# The player can combine every action with the opportunity
		# to manipulate towards a certain signal
		g1 = g2 = matrix(NA,na[1],na[2])
		
		for (y1 in 0:m$ny) {
			cost1 = cost * (y1>0)
			for (y2 in 0:m$ny) {
				cost2 = cost * (y2>0)
				rows = (y1*ona[1]+1):(y1*ona[1]+ona[1])
				cols = (y2*ona[2]+1):(y2*ona[2]+ona[2])
	
				g1[rows,cols]=og1-cost1
				g2[rows,cols]=og2-cost2
			}
		}
		
		rownames(g1) = rownames(g2) = lab.a1
		colnames(g1) = colnames(g2) = lab.a2

		phi.mat = get.dummy.phi.mat(ny,prod(na))
		rownames(phi.mat) = m$lab.y
		mm = init.game(g1=g1,g2=g2,phi.mat=phi.mat)
		phi.mat = mm$phi.mat
		
		
		###################################################################################
		# Creating mapping of action profiles of manipulation game and original game
		###################################################################################
		
		a1.map = cbind(1:na[1],make.grid.matrix(x=list(0:m$ny,1:ona[1])))
		a2.map = cbind(1:na[2],make.grid.matrix(x=list(0:m$ny,1:ona[2])))
		colnames(a1.map)=colnames(a2.map)= c("ai","my","oai")

		
		ai.mat = get.ai.mat(mm)
		i = 1
		map = matrix(NA,mm$nA,8)
		colnames(map)=c("a","a1","a2","oa","oa1","oa2","my1","my2")
		map[,"a"] = 1:mm$nA; 
		map[,c("a1","a2")]=get.ai.mat(mm);
		map[,"oa1"] = a1.map[map[,"a1"],"oai"]
		map[,"oa2"] = a2.map[map[,"a2"],"oai"]
		map[,"oa"]  = get.a.from.ai.mat(m,map[,c("oa1","oa2")])

		map[,"my1"]  = a1.map[map[,"a1"],"my"]
		map[,"my2"]  = a2.map[map[,"a2"],"my"]

		####################################################################################
		# Construct new phi.mat
		####################################################################################
		
		phi.mat=m$phi.mat[,map[,"oa"]]
		man.prob1 = man.prob2 = matrix(0,NROW(phi.mat),NCOL(phi.mat))
		rows = map[,"my1"] > 0
		man.prob1[map[rows,c("my1","a")]]=1
		man.prob2[map[rows,c("my2","a")]]=1
		colnames(man.prob1) = colnames(man.prob2)= colnames(phi.mat)
		man.prob2
		colnames(phi.mat)=mm$lab.a; rownames(phi.mat)= mm$lab.y
		
		sum.man = apply(man.prob1,2,sum)+apply(man.prob2,2,sum)
		sum.man = matrix(sum.man,NROW(phi.mat),NCOL(phi.mat),byrow=TRUE)
								
		phi.mat = (1 - pi.man*sum.man)*phi.mat + man.prob1*pi.man + man.prob2*pi.man
		colnames(phi.mat)=mm$lab.a; rownames(phi.mat)= mm$lab.y
		phi.mat
		mm$phi.mat = phi.mat
		
		return(mm)
}

