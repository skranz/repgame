# R code of the examples analysed in 
# "Interactively Solving Repeated Games: A Toolbox"
#
# Author: Sebastian Kranz, Ulm University
# Version: 0.3
# Date: 03.11.2014
#
# Examples of Part 2 Sections 7-10



##################################################################
# Section 7: 
# Noisy Prisoners' Dilemma Game
##################################################################


library(repgame)

init.noisy.pd.game = function(d,s,lambda,mu,alpha,beta,psi) {
  store.objects("init.noisy.pd.game")
  # restore.objects("init.noisy.pd.game")

  # Payoff matrix of player 1 (that of player 2 is symmetric)
  g1 = matrix(c(
    1,    -s,
    1+d,   0),2,2,byrow=TRUE)
  lab.ai = c("C","D")
      
  # Create phi.mat, which stores the signal distributions
  prob.yC = c(1-2*alpha-lambda,         1-2*alpha-lambda-mu-beta,
              1-2*alpha-lambda-mu-beta, 1-psi)
                           
  phi.mat = matrix(c(
    #CC       CD          DC         DD 
    prob.yC,                                #yC
    lambda,   lambda+mu,  lambda+mu , psi,  #yD
    alpha,    alpha,      alpha+beta, 0,    #y1
    alpha,    alpha+beta, alpha     , 0     #y2
  ), 4,4,byrow=TRUE)
  lab.y = c("yC","yD","y1","y2")
  m=init.game(g1=g1,phi.mat=phi.mat, symmetric=TRUE, 
              lab.ai = lab.ai,lab.y=lab.y, name="Noisy PD")
  m
}

d=1;s=1.5;lambda=0.1;mu=0.2;alpha=0.1;beta=0.2;psi = 1
m = init.noisy.pd.game(d,s,lambda,mu,alpha,beta,psi)
m = solve.game(m)
plot(m,xvar="L",yvar=c("Ue","V"),lwd=2)
plot(m,xvar="delta",delta.seq=seq(0.5,1,by=0.001),
     yvar=c("Ue","V"),lwd=2)

##################################################################
# Section 8: 
# Green-Porter Style Models
##################################################################


     
## Section 8.2 ###################################################
     
     
#  Require  global  parameters  to  be  set
#  D0.fun,  p.max,  F.s  =  F.s,  cost.fun  =  cost.fun,q.max  =  q.max
#  Either  nq=nq  or  q.seq=NULL
init.payoffs.green.porter  =  function(make.q.seq=TRUE)  {
  restore.point("init.payoffs.green.porter")
  
  #  Generate  the  distribution  function  for  prices  F.p
  #  The  assignment  <<-  stores  F.p  in  the  global  environment
  F.p  <<-  function(p,Q)  {
    x  =  F.s(log(Q  /  D0.fun(p)))
    x[p  >=  p.max]  ==  0
    x
  }
  #  Action  space  of  each  firm  (we  assume  firms  are  symmetric)
  if  (make.q.seq)  {
    q.seq  <<-  seq(0,q.max,length=nq)
  }
  #  Generate  function  that  calculates  expected  payoffs  for
  #  each  action  profile
  #  The  function  performs  a  very  fine  discretization  of
  #  the  distribution  F.p  to  approximate  expected  payoffs
  g.fun  =  function(qm)  {
    store.objects("g.fun")
    #restore.objects("g.fun")
    Q  =  rowSums(qm)
    g  =  matrix(0,NROW(qm),2)
    #  Calculate  expected  profits  for  every  action  profile
    #  The  loop  takes  some  time,  since  calculation  of  an  expected  price
    #  uses  numerical  integration  for  every  level  of  Q
    for  (r  in  1:NROW(qm))  {
      Ep  =  calc.mean.from.F.fun(F.p,Q=Q[r],x.min=0,x.max=p.max)
      g[r,]  =  Ep  *  qm[r,]  -  cost.fun(qm[r,])
    }
    #  Manually  correct  the  case  that  no  firm  produces  anything
    if  (Q[1]==0)  {
      g[1,]  =  c(0,0)
    }
    g
  }
  m  =  init.game(g.fun=g.fun,  action.val  =  q.seq,  symmetric=TRUE,
  name="Green-Porter  with  Optimal  Penal  Codes")
  return(m)
}

init.signals.green.porter  =  function(m=m)  {
  qm  =  m$action.val.mat
  Q=rowSums(qm)
  y.val  <<-  seq(0,p.max,length=ny)
  #  mapply  and  apply  work  similar  than  Vecorize
  #  Payoffs  at
  F.y  =  mapply(F.p,Q=Q,MoreArgs=list(p=y.val))
  phi.mat  =  apply(F.y,2,discretize.given.F.vec)
  phi.mat[,Q==0]  =  c(rep(0,ny-1),1)
  colSums(phi.mat)
  m  =  set.phi.mat(m,phi.mat,lab.y  =  paste("y",y.val,sep=""))
  m$y.val  =  y.val
  return(m)
}

init.signals.green.porter  =  function(m=m)  {
  qm  =  m$action.val.mat
  Q=rowSums(qm)
  y.val  <<-  seq(0,p.max,length=ny)
  

  #  Calculate  the  F.p  at  every  signal  y
  F.y  =  mapply(F.p,Q=Q,MoreArgs=list(p=y.val))
  
  #  Assign  probability  weight  according  to  the
  #  discretization  procedure  explained  in  the  tutorial
  phi.mat  =  apply(F.y,2,discretize.given.F.vec)
  
  #  We  have  to  manually  correct  the  price  distribution  if
  #  no  firm  produces  anything
  phi.mat[,Q==0]  =  c(rep(0,ny-1),1)
  
  #  Initialize  the  signal  distribution  for  the  model  m
  m  =  set.phi.mat(m,phi.mat,lab.y  =  paste("y",y.val,sep=""))
  m$y.val  =  y.val
  return(m)
}
#  Initialize  parameters

#  Deterministic  part  of  demand  function
alpha  =  100;  beta  =  1;
D0.fun  =  function(P)  {alpha-beta*P}
p.max  =  alpha  /  beta
q.max  =  alpha


#  Distribution  of  market  size
mean.s  =  0;  sigma.s  =  0.2
F.s  =  function(x)  {pnorm(x,mean.s,sigma.s)}

#  Cost  function
MC  =  10
cost.fun  =  function(q)  {MC*q}

#  Number  of  actions  per  firm  and  number  of  signals
nq  =  51;
ny  =  21;

nq  =  21;
ny  =  11;


#  Initialize  model  with  perfect  monitoring
m.pm  =  init.payoffs.green.porter()
m.pm$name  =  "Cournot  with  Payments  &  Optimal  Penal  Codes"

#  Transform  into  a  model  with  imperfect  monitoring
m  =  init.signals.green.porter(m.pm)
m$name  =  "Green-Porter  with  Payments  &  Optimal  Penal  Codes"


## Section 8.3 ###################################################
     

#  Solve  the  models of perfect and imperfect monitoring
m.pm = solve.game(m.pm,  ignore.ae  =  m$action.val[,1]>m$action.val[,2])
m = solve.game(m,  ignore.ae  =  m$action.val[,1]>m$action.val[,2])


# Compare models graphically
plot.compare.models(m.pm,m,xvar="delta",yvar="Ue",
          m1.name = "Cournot",m2.name="Green-Porter Style",
          legend.pos ="bottomright", identify="Ue")



# Analyse strategy profiles
ret  =  levelplot.rep(m=m,z=m$G,main="G  and  opt.  action  plan",
          focus=2,cuts=100,  xlab="q1",ylab="q2",  identify=NULL,
          col.nash="purple",col.br1  =  "green",  col.br2  =  "greenyellow",
          col.ae  =  "blue",  col.a1  =  "red")
 

## Section 8.4 Optimal equilibrium action profiles ############################
          
ret = levelplot.rep(m=m,z=m$G,main="G and opt. ae",focus=2,cuts=100,
                    xlab="q1",ylab="q2", identify="Ue",
                    xlim=c(-1,35),ylim=c(20,56), zlim=c(500,2000),
                    col.nash="purple",col.br1 = NA, col.br2 = NA,
                    col.ae = "blue", col.a1 = NA)
                    

                    
                    
         
# Solving the game with symmetric action profiles
ai.mat  =  get.ai.mat(m)
ignore.ae  =  !(ai.mat[,1]==ai.mat[,2]  |
ai.mat[,1]  ==  (ai.mat[,2]-1)  )
m.sym.ae  =  solve.game(m,ignore.ae=ignore.ae)


plot.compare.models(m,m.sym.ae,xvar="L",yvar="Ue",
          m1.name = "unrestricted", m2.name = "sym. ae",
          legend.pos  =  "bottomright",identify="Ue")

          
# Section 8.5 ####################################################

# Level plots for the punishment state

ret = levelplot.rep(m=m,z=m$c[,1],main="c1 and opt.a1",
  focus=0,cuts=100,
  xlab="q1",ylab="q2", identify="v1",
  xlim=c(-1,35),ylim=c(25,100),
  col.nash="purple",col.br1 = "green", col.br2 = "greenyellow",
  col.ae = NA, col.a1 = "red")
  

# Show liquidity requirements

ret = levelplot.rep(m=m,z=m$L,main="L",focus=-3,cuts=100,
  xlab="q1",ylab="q2", identify="v1",
  xlim=c(-1,35),ylim=c(25,100),
  col.nash="purple",col.br1 = "green", col.br2 = "greenyellow",
  col.ae = NA, col.a1 = NA)
  

# Solve game by assuming that player 1 plays a best reply in his punishment  
ignore.a1 = (m$g[,1] != m$c[,1])                                        
m.a1.br1 = solve.game(m,ignore.a1=ignore.a1)
plot.compare.models(m,m.a1.br1,xvar="L",yvar="v1", x.max=4000,
                    m1.name="m", m2.name="m1.a1.br1")
                    
       
#################################################################################
# Section 9
#################################################################################

# Analyse Green-Porter with auditing  and compare to case without audting
         
# Solving the game with symmetric action profiles


#  Initialize  model  with  perfect  monitoring
m.pm  =  init.payoffs.green.porter()
m.imp  =  init.signals.green.porter(m.pm)
m.10 = set.audit.prob(m.imp,0.1)


m.imp$name  =  "Green-Porter (no audits)"
m.pm$name = "Green-Porter (permanent audits)"
m.10$name = "Green-Porter (10% audits)"
           
ai.mat  =  get.ai.mat(m)
ignore.ae  =  !(ai.mat[,1]==ai.mat[,2]  |
ai.mat[,1]  ==  (ai.mat[,2]-1))
m.10   = solve.game(m.10 ,ignore.ae=ignore.ae)
m.pm   = solve.game(m.pm ,ignore.ae=ignore.ae)
m.imp  = solve.game(m.imp,ignore.ae=ignore.ae)

plot.compare.models(m.10,m.imp,xvar="delta",yvar="Ue",
                    m1.name = "10% audit", m2.name = "no audits",
                    legend.pos ="bottomright")

plot.compare.models(m.10,m.pm,xvar="delta",yvar="Ue",
                    m1.name = "10% audit", m2.name = "100% audits",
                    legend.pos ="bottomright")

                    
  