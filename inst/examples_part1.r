# R code of the examples analysed in 
# "Interactively Solving Repeated Games: A Toolbox"
#
# Author: Sebastian Kranz, Ulm University
# Version: 0.3
# Date: 03.11.2014
#

# First install the packackes. See the tutorial


##################################################################
# Section 3: 
# Abreu's simple Cournot Game
##################################################################

# Define Payoff Matrices for player 1 and 2
              #L #M  #H
g1 = rbind(c( 10, 3, 0),    # L
           c( 15, 7, -4),   # M
           c(  7, 5,-15))   # H
g2 = t(g1)  # Player 2's payoff matrix is simply the transpose of that of player 1

# Load package

library(repgame)
# Initialize the model of the repeated game
m = init.game(g1=g1,g2=g2,lab.ai=c("L","M","H"),
              name="Simple Cournot Game",symmetric=FALSE)
# Solve the model 
m = solve.game(m, keep.only.opt.rows=!TRUE)

# Show matrix with critical delta, optimal action structures and payoffs
m$opt.mat

# Plot optimal payoffs and show matrix that characterizes optimal strategy profiles and the cutoffs of the discount factor
plot(m,xvar="delta",yvar=c("Ue","V"))


# Try the same plot with interactive mode. You can click with the left mouse
# button on the plot
plot(m,xvar="delta",yvar=c("Ue","V"),identify=TRUE)


# Solve the model with grim-trigger strategies
m.gt = set.to.grim.trigger(m)
m.gt = solve.game(m.gt)
m.gt$opt.mat

plot(m.gt)

# Plot comparision of optimal equilibria and grim trigger equilibria
plot.compare.models(m,m.gt,xvar="delta",yvar="Ue",legend.pos="topleft",
                    m1.name = "opt", m2.name="grim",identify=TRUE)

                    
##################################################################
# Section 4: Public Goods Games
##################################################################

### Section 4.1 ##################################################
                    
# Function that initializes, solves and plots a public goods game
public.goods.game = function(X, k1,k2=k1) {
  # Two lines that are useful for debugging (see hints below)
  restore.point("public.goods.game")
   
  # Initialize matrices of players contributions
  x1m = matrix(X,NROW(X),NROW(X),byrow=FALSE) 
  x2m = matrix(X,NROW(X),NROW(X),byrow=TRUE)                      

  # Calculate payoff matrices for each player
  g1 = (x1m+x2m)/2 - k1*x1m         
  g2 = (x1m+x2m)/2 - k2*x2m
  
  # Give the game a name
  name = paste("Public Goods Game k1=", k1, " k2=",k2,sep="")

  # Initialize the game	    
  m = init.game(n=2,g1=g1,g2=g2,name=name, lab.ai=X)
  
  # Solve and plot the model 
  m = solve.game(m, keep.only.opt.rows=TRUE)
  plot(m)
  return(m)
}

# Call the function
m = public.goods.game(X=0:10,k1=0.6,k2=0.75)
# Look at the optimal solution
m$opt.mat



get.mat.of.delta(m,delta=c(0.2,0.8))

### Section 4.2 Comparative Statics  ################################################

# 4.2.1                    
pg.games = Vectorize(public.goods.game, vectorize.args = c("k1","k2"),SIMPLIFY=FALSE)
k.seq = seq(0.5,1,by = 0.01)
m.list = pg.games(X=0:10,k1=k.seq,k2=k.seq)

mat = levelplot.payoff.compstat(m.list, par = k.seq ,
        xvar = "k", yvar = "delta", payoff.var="Ue",
        delta = seq(0,1,by = 0.01),col.scheme = "grey")

# 4.2.2
 
# Generate grid of k1 and k2 combinations between 0.5 and 1
k.seq = seq(0.5,1,by = 0.02)
  
# (make.grid.matrix is an own function similar to expand.grid)
k.grid = make.grid.matrix(x=k.seq,n=2)
colnames(k.grid)=c("k1","k2")
  
# Solve the model for all combinations of k1 and k2
m.list = pg.games(X=0:10,k1=k.grid[,1],k2=k.grid[,2])
  
# Show level plot of payoffs for delta=0.5
delta = 0.5
levelplot.payoff.compstat(m.list, par = k.grid , 
  xvar = "k1", yvar = "k2", payoff.var="Ue",delta = delta, 
  col.scheme = "grey")

# Generate a sequence of levelplots for different discount factors
for (delta in seq(0,1,by=0.1)) {
  levelplot.payoff.compstat(m.list, par = k.grid, 
    xvar = "k1", yvar = "k2", payoff.var="Ue",delta = delta, 
    col.scheme = "grey") 
}



### Section 4.3 n-player public goods games  ####################################


# Creates a n-player public goods game
# n is the number of players. 
# X is the set of different contribution levels
# k is a n x 1 vector with production costs for every player
pg.game = function(n,X,k=rep(((1+1/n)/2),n)) {
  # A function that returns a payoff matrix
  # given a action matrix x.mat, of which
  # every row corresponds to one action profile
  g.fun = function(x.mat) {
    g = matrix(0,NROW(x.mat),n)
    xsum = rowSums(x.mat)
    for (i in 1:n) {
      g[,i] = xsum / n - k[i]*x.mat[,i]
    }
    g
  }
  name=paste(n,"Player Public Goods Game")
  m = init.game(n=n,g.fun=g.fun,action.val = X,
	             name=name, lab.ai=round(X,2))
  m=solve.game(m)
  plot(m)
  m
}
# Solve a 3 player public goods game with
# asymmetric production costs
m = pg.game(n=3,X=0:10,k=c(0.4,0.6,0.8))

# Let us have a look at the payoffs as function of L
plot(m,xvar="L")


             
##################################################################
# Section 5: Cournot Games
##################################################################

# Section 5.1 ##############################################
init.cournot = function(n=2,q.seq, A, B, MC) {
  P.fun = function(Q) {A-B*Q}
  cost.fun = function(q) {MC*q}  

  g.fun = function(qm) {
    Qm = rowSums(qm)
    Pm = P.fun(Qm)
    Pm[Pm<0]=0
    g = matrix(NA,NROW(qm),n)
    for (i in 1:n) {
      g[,i] = Pm*qm[,i]-cost.fun(qm[,i])
    }
    g
  }
  name=paste(n,"Cournot Game (g.fun)")
  m = init.game(n=2,g.fun=g.fun,action.val = q.seq,
                name=name, lab.ai=round(q.seq,2))
  return(m)
}
# Parameters of the model
n=2; A=100; B=1; MC=10;
q.seq = seq(0,100,by=1); # Grid of possible quantities
m = init.cournot(n=n,q.seq=q.seq,A=A,B=B,MC=MC)
m = solve.game(m)
plot(m,legend.pos = "right")


# Section 5.2 ###########################################################

pmx.init.cournot = function(n=2, q.range, A, B, MC) {

  # First part
  store.objects("pmx.init.cournot")
  #restore.objects("pmx.init.cournot")

  # Stage game payoffs
  gx.fun = function(q) {
    Q = rowSums(q)
    P = A-B*Q
    P[P<0]=0
    g = matrix(NA,NROW(q),n)
    for (i in 1:n) {
      g[,i] = (P-MC)*q[,i]
    }
    g
  }

  # Stage game cheating payoffs
  cx.fun = function(q) {
    Q   = rowSums(q)
    c.mat = matrix(NA,NROW(q),n)
    for (i in 1:n) {
      Q_i = Q-q[,i]
      c.mat[,i] = (A-MC-B*(Q-q[,i]))^2 / (4*B)
      c.mat[A-MC-B*(Q-q[,i])<0,i] = 0
    }
    c.mat
  }
  
  # Second part
  name=paste("Cournot Game (",n," firms)",sep="")
  m = pmx.init.game(n=n,nx=n,x.range=q.range,symmetric=TRUE,
        gx.fun=gx.fun,cx.fun=cx.fun, 
        name=name, names.x=paste("q",1:n,sep=""))
  return(m)
} 

# Parameters of the model
n=4; A=100; B=1; MC=10;
q.max = A / (n-1)   # maximal output a firm can choose
m = pmx.init.cournot(n=n,q.range=c(0,q.max),A=A,B=B,MC=MC)

cnt = list(method="grid", step.size.start=5, step.size.end = 0.5,
          num.refinements = 3)
m = solve.game(m,cnt=cnt)
m.grid = m

# Solve with adaptive grid refinement
cnt = list(method="grid",
           num.step.start = 8, num.step.end = 128 ,
           step.size.factor = 4, grid.size = 100000)
m.grid = solve.game(m,cnt=cnt)


# Solve with random sampling
cnt = list(method="random", step.size=0.1,
           size.draw = 1000,num.draws = 200,
           local.abs.width=4, prob.draw.global = 0.3)
m.ran = solve.game(m,cnt=cnt)

# Compare the two models graphically
plot.compare.models(m.grid,m.ran, xvar="delta", yvar="Ue",
                    m1.name="grid", m2.name="random", 
                    identify = TRUE, legend.pos="right")

# Section 5.3          ###########################################

pmx.init.cournot.sym = function(n=2, q.range, A, B, MC) {

  # First part
  store.objects("pmx.init.cournot.sym")
  #restore.objects("pmx.init.cournot.sym")

  # Stage game payoffs
  gx.fun = function(q) {
    Q = rowSums(q)
    P = A-B*Q
    P[P<0]=0
    g = matrix(NA,NROW(q),n)
    for (i in 1:n) {
      g[,i] = (P-MC)*q[,i]
    }
    g
  }

  # Stage game cheating payoffs
  cx.fun = function(q) {
    Q   = rowSums(q)
    c.mat = matrix(NA,NROW(q),n)
    for (i in 1:n) {
      Q_i = Q-q[,i]
      c.mat[,i] = (A-MC-B*(Q-q[,i]))^2 / (4*B)
      c.mat[A-MC-B*(Q-q[,i])<0,i] = 0
    }
    c.mat
  }
  

  # Second part: Code that specifies symmetry constraints

  # In equilibrium state only player 1 chooses freely his
  # action other firms are assumed to choose symmetric actions
  x.free.e = 1 
  link.fun.e = function(xf.mat,k=NULL) {
    # xf.mat will only have one column, which 
    # contains the actions of firm 1
    # return a matrix with n cols that are all identical to xf.mat
    return(matrix(as.vector(xf.mat),NROW(xf.mat),n))
  }
  # In the punishment states only the punished player and one other player
  # (here either player 1 or 2) can freely decide on their actions
  x.free.i = lapply(1:n, function(i) c(1,i))
  x.free.i[[1]] = c(1,2)
  link.fun.i = function(xf.mat,k) {
    if (k == 1) {
      mat = matrix(as.vector(xf.mat[,2]),NROW(xf.mat),n)
      mat[,1] = xf.mat[,1] 
    } else {
      mat = matrix(as.vector(xf.mat[,1]),NROW(xf.mat),n)
      mat[,k] = xf.mat[,k] 
    }
    return(mat)
  }
  m = pmx.init.game(n=n,nx=n,gx.fun=gx.fun,cx.fun=cx.fun,
              x.range=q.range, symmetric=TRUE,
              x.free.e = x.free.e, link.fun.e=link.fun.e, 
              x.free.i = x.free.i,link.fun.i = link.fun.i,
              name=paste("Cournot Game (",n," firms)",sep=""),
              names.x=paste("q",1:n,sep=""))
  return(m)
}

# Parameters of the model
n=100; A=100; B=1; MC=10;
m = pmx.init.cournot(n=n,q.range=c(0,A / (n-1)),A=A,B=B,MC=MC)



##################################################################
# Section 6: Hotelling Games
##################################################################


#6.2 The model with a wrong best-reply function  
init.hotelling = function(ms,symmetric=FALSE) {
  colnames(ms) = c("market","firm1","firm2","MC1","MC2","tau","size","w")
  n = max(ms[,c("firm1","firm2")])
  
	restore.point("init.multimarket.hotelling")
	
  # Stage game payoffs
  gx.fun = function(p) {
    store.objects("gx.fun")
    #restore.objects("gx.fun")
    g.mat = matrix(0,NROW(p),n)    
    for (ma in 1:NROW(ms)) {
      tau = ms[ma,"tau"];size = ms[ma,"size"]; w = ms[ma,"w"];
      p1 = p[,ma*2-1];p2 = p[,ma*2];
      q1 = pmax(0,pmin((w-p1)/tau,pmin(1,(1/2+(p2-p1)/(2*tau)))))*size
      q2 = pmax(0,pmin((w-p2)/tau,pmin(1,(1/2+(p1-p2)/(2*tau)))))*size
      g.mat[,1] = g.mat[,1]+q1*(p1-ms[ma,"MC1"])
      g.mat[,2] = g.mat[,2]+q2*(p2-ms[ma,"MC2"])
    }
    g.mat
  }
  
  # Stage game cheating payoffs. Assuming FOC specfies best reply prices
  cx.fun = function(p) {
    store.objects("cx.fun")
    #restore.objects("cx.fun");restore.objects("init.multimarket.hotelling")
    c.mat = matrix(0,NROW(p),n)
    ma = 1
    while(ma <= NROW(ms)) {
      tau = ms[ma,"tau"];size = ms[ma,"size"];w = ms[ma,"w"];
      MC1 = ms[ma,"MC1"]; MC2 = ms[ma,"MC2"];
      p1 = p[,ma*2-1];p2 = p[,ma*2];
      p1.br = (1/2)*(p2+MC1+tau)
      p2.br = (1/2)*(p1+MC2+tau)
      q1 = pmax(0,pmin((w-p1.br)/tau,pmin(1,(1/2+(p2-p1.br)/(2*tau)))))*size
      q2 = pmax(0,pmin((w-p2.br)/tau,pmin(1,(1/2+(p1-p2.br)/(2*tau)))))*size
     
      c.mat[,1] = c.mat[,1]+q1*(p1.br-MC1)
      c.mat[,2] = c.mat[,2]+q2*(p2.br-MC2)
      ma = ma+1
    }
    
    c.mat
  }
  name=paste("Hotelling (markets: ",NROW(ms),")", sep="")
  nx = NROW(ms)*2
  name.grid = make.grid.matrix(x=list(1:NROW(ms),1:2))
  names.x = paste("m",name.grid[,1],".p",name.grid[,2],sep="")
  x.range = cbind(0,rep(ms[,"w"],each=2))
  m = pmx.init.game(n=n,nx=nx,x.range=x.range,symmetric=symmetric,
        gx.fun=gx.fun,cx.fun=cx.fun, 
        name=name, names.x=names.x)
  return(m)
} 
ms = matrix(c(
 #market  firm1, firm2,   MC1,  MC2, tau, size, w,
       1,    1,     2,      10,  10, 100, 100,   200
     ),ncol=8,byrow=TRUE)
m = init.hotelling(ms=ms,symmetric=TRUE)
cnt = list(method="grid", step.size.start=10,step.size.end = 0.5,
          num.refinements = 1)
m = solve.game(m,cnt=cnt)
plot(m,legend.pos="right")
m.unsure = m


mat = check.cx.fun(m,num=30,method="grid",num.grid.steps=10000)
round(mat[,"diff"],3)

mat = check.cx.fun(m,x.sample="opt.mat",
                   method="grid",num.grid.steps=10000)
round(mat[,"diff"],3)

# Section 6.3 specifying cheating payoffs from numerical optimization
m = init.hotelling(ms=ms,symmetric=TRUE)
m$cx.fun = make.numerical.cx.fun(m,method="grid",num.grid.steps=1000)
m = solve.game(m,cnt=cnt)
plot(m,legend.pos="right")
m.num = m

plot.compare.models(m.unsure,m.num)

# Section 6.4 correct best-reply functions ######################################################

 
init.hotelling = function(ms,symmetric=FALSE) {
  colnames(ms) = c("market","firm1","firm2","MC1","MC2","tau","size","w")
  n = max(ms[,c("firm1","firm2")])
  
  store.objects()
	#restore.objects("init.multimarket.hotelling")
	
  # Stage game payoffs
  gx.fun = function(p) {
    store.objects("gx.fun")
    #restore.objects("gx.fun")
    g.mat = matrix(0,NROW(p),n)    
    for (ma in 1:NROW(ms)) {
      tau = ms[ma,"tau"];size = ms[ma,"size"]; w = ms[ma,"w"];
      p1 = p[,ma*2-1];p2 = p[,ma*2];
      q1 = pmax(0,pmin((w-p1)/tau,pmin(1,(1/2+(p2-p1)/(2*tau)))))*size
      q2 = pmax(0,pmin((w-p2)/tau,pmin(1,(1/2+(p1-p2)/(2*tau)))))*size
      g.mat[,1] = g.mat[,1]+q1*(p1-ms[ma,"MC1"])
      g.mat[,2] = g.mat[,2]+q2*(p2-ms[ma,"MC2"])
    }
    g.mat
  }
  
  # Stage game cheating payoffs. Assuming FOC specfies best reply prices
  cx.fun = function(p) {
    store.objects("cx.fun")
    #restore.objects("cx.fun");restore.objects("init.multimarket.hotelling")
    c.mat = matrix(0,NROW(p),n)
    ma = 1
    while(ma <= NROW(ms)) {
      tau = ms[ma,"tau"];size = ms[ma,"size"];w = ms[ma,"w"];
      
      for (i in 1:2) {
        j = 3-i
        # Not very elegant but who cares...
        if (i == 1) {
          MC.i = ms[ma,"MC1"]; MC.j = ms[ma,"MC2"];
          p.i = p[,ma*2-1];p.j = p[,ma*2];
        } else {
          MC.j = ms[ma,"MC1"]; MC.i = ms[ma,"MC2"];
          p.j = p[,ma*2-1];p.i = p[,ma*2];
        }
        
        # Try out the different best replies

                
        # kink
        p.i.br = 2*w-tau-p.j
        qi = pmax(0,pmin((w-p.i.br)/tau,pmin(1,(1/2+(p.j-p.i.br)/(2*tau)))))*size
        pi.k = qi*(p.i.br-MC.i)
        
        # interior
        p.i.br = ((tau+MC.i+p.j)/2)
        qi = pmax(0,pmin((w-p.i.br)/tau,pmin(1,(1/2+(p.j-p.i.br)/(2*tau)))))*size
        pi.c = qi*(p.i.br-MC.i)
        
        # steal the whole market from the other player
        p.i.br = p.j-tau
        qi = pmax(0,pmin((w-p.i.br)/tau,pmin(1,(1/2+(p.j-p.i.br)/(2*tau)))))*size
        pi.s = qi*(p.i.br-MC.i)
 
        # Set monopoly price in those rows where pj is very high
        pi.m = rep(0,NROW(p))
        rows = (p.j >= (3/2)*w-tau-(1/2)*MC.i)
        p.i.br = (w+MC.i)/2
        qi = pmax(0,pmin((w-p.i.br)/tau,pmin(1,(1/2+(p.j-p.i.br)/(2*tau)))))*size
        pi.m[rows] = (qi*(p.i.br-MC.i))[rows]
               
        # Take that best-reply that yields highest profits
        c.mat[,i] = c.mat[,i]+pmax(pi.k,pmax(pi.c,pmax(pi.s,pi.m)))
      }
      ma = ma +1
    }
    return(c.mat)
  }
  name=paste("Hotelling (markets: ",NROW(ms),")", sep="")
  nx = NROW(ms)*2
  name.grid = make.grid.matrix(x=list(1:NROW(ms),1:2))
  names.x = paste("m",name.grid[,1],".p",name.grid[,2],sep="")
  x.range = cbind(0,rep(ms[,"w"],each=2))
  m = pmx.init.game(n=n,nx=nx,x.range=x.range,symmetric=symmetric,
        gx.fun=gx.fun,cx.fun=cx.fun, 
        name=name, names.x=names.x)
  return(m)
} 


ms = matrix(c(
 #market  firm1, firm2,   MC1,  MC2, tau, size, w,
       1,    1,     2,      10,  10, 100, 100,   200
     ),ncol=8,byrow=TRUE)
m = init.hotelling(ms=ms,symmetric=TRUE)
cnt = list(method="grid", step.size.start=10,step.size.end = 0.5,
          num.refinements = 1)
m = solve.game(m,cnt=cnt)
plot(m,legend.pos="right")


  

mat = check.cx.fun(m,num=100,method="grid",num.grid.steps=10000)
round(mat[,"diff"],3)

plot.compare.models(m,m.unsure)


# Perform comparative statics in tau
m.list = list()
tau.seq = seq(1,120,by=5)
cnt = list(method="grid", step.size.start=10,step.size.end = 0.5,
          num.refinements = 1)
for (k in 1:NROW(tau.seq)) {
  tau = tau.seq[k]
ms = matrix(c(
 #market  firm1, firm2,   MC1,  MC2, tau, size, w,
       1,    1,     2,      10,  10, tau, 100,   200
     ),ncol=8,byrow=TRUE)
  m = init.hotelling(ms=ms,symmetric=TRUE)
  m$name = paste("tau: ", tau,sep="")
  m = solve.game(m,cnt=cnt)
  plot(m,legend.pos="right")
  m.list[[k]] = m
}
mat = levelplot.payoff.compstat(m.list, par = tau.seq ,
        xvar = "tau", yvar = "delta", payoff.var="Ue",
        delta = seq(0,0.6,by = 0.01),col.scheme = "grey", identify=TRUE)

# Perform comparative statics in tau and MC1
m.list = list()
tau.seq = c(10,25,50,100)
MC1.seq = c(0,10,25,50,100)
par.grid = make.grid.matrix(x=list(tau.seq,MC1.seq))
colnames(par.grid)=c("tau","MC1")

for (k in 1:NROW(par.grid)) {
  tau = par.grid[k,1]
  MC1 = par.grid[k,2]
ms = matrix(c(
 #market  firm1, firm2,   MC1,  MC2, tau, size, w,
       1,    1,     2,      MC1,  0, tau, 100,   200
     ),ncol=8,byrow=TRUE)
  m = init.hotelling(ms=ms,symmetric=TRUE)
  m$name = paste("tau: ", tau, " MC1: ", MC1,sep="")
  m = solve.game(m,cnt=cnt)
  plot(m,legend.pos="right")
  m.list[[k]] = m
}

for (delta in seq(0,1,by = 0.1)) {
  mat = levelplot.payoff.compstat(m.list, par = par.grid ,
        xvar = "tau", yvar = "MC1", payoff.var="Ue",
        delta = delta,col.scheme = "grey", identify=NULL)
}


# Section 6.5 Multimarket ######################################################

 
# First let us compare two single markets
# One with high transportation cost the other with low

cnt = list(method="grid", step.size.start=10,step.size.end = 1, num.refinements = 1)

ms = matrix(c(
 #market  firm1, firm2,   MC1,  MC2, tau, size, w,
       1,    1,     2,      10,  10, 100, 100,   400
     ),ncol=8,byrow=TRUE)
m = init.hotelling(ms=ms,symmetric=TRUE)
m.high = solve.game(m,cnt=cnt)



ms = matrix(c(
 #market  firm1, firm2,   MC1,  MC2, tau, size, w,
       1,    1,     2,      10,  10, 10,  100,   400
     ),ncol=8,byrow=TRUE)
m = init.hotelling(ms=ms,symmetric=TRUE)
m.low = solve.game(m,cnt=cnt)

plot.compare.models(m.high,m.low,xvar="delta",yvar="Ue",
                    m1.name="high tau",m2.name="low tau",
                    legend.pos = "bottomright")
                    


# Parameters of the model
ms = matrix(c(
 #market  firm1, firm2,   MC1, MC2, tau, size, w,
       1,    1,     2,      10,  10, 100, 100,   400,
       1,    1,     2,     10,   10,  10, 100,   400
     ),ncol=8,byrow=TRUE)

m = init.hotelling(ms=ms,symmetric=TRUE)

cnt = list(method="grid", step.size.start=5,step.size.end = 2, num.refinements = 1,
           use.random.start.points = TRUE, num.random.start.points=10000)           
m = solve.game(m,cnt=cnt)

# As a check, try to refine the solution with random sampling
cnt = list(method="random", step.size=1,
           size.draw = 100,num.draws = 50,
           local.abs.width=4, prob.draw.global = 1)
m = refine.solution(m,cnt)


m.mult = m
plot(m.mult)

mat = check.cx.fun(m.mult,i.of.x = c(1,2,1,2), use.i = c(1,2),
                   method="grid",num.grid.steps=100)


# Compare the two models graphically

# Generate matrices for same delta.sequence
delta = seq(0,1,by=0.01)
mat1 = get.mat.of.delta(m.high,delta = delta,add.crit.delta=FALSE)
mat2 = get.mat.of.delta(m.low,delta = delta,add.crit.delta=FALSE)
mat.b = get.mat.of.delta(m.mult,delta = delta,add.crit.delta=FALSE)
mat.s  = mat1
mat.s[,"Ue"]= mat1[,"Ue"]+mat2[,"Ue"]
# Draw the different matrices
ylim = range(c(mat.s[,"Ue"],mat.b[,"Ue"]))
plot(mat.b[,"delta"],mat.b[,"Ue"],col="blue",type="l",
     main="Single vs Multi Market",xlab="delta",ylab="Ue",ylim=ylim)
lines(mat.s[,"delta"],mat.s[,"Ue"],col="red",lty=1)
lines(mat.b[,"delta"],mat.b[,"Ue"],col="blue",lty=2)
lines(mat1[,"delta"],mat1[,"Ue"],col="green",lty=2)
lines(mat2[,"delta"],mat2[,"Ue"],col="green",lty=2)

