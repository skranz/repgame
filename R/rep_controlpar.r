# R code for the package repgames
#
# Author: Sebastian Kranz, University of Bonn
#         skranz@uni-bonn.de
#          
# 
# "Interactively Solving Repeated Games: A Toolbox"



# Has to be called once before the model can be used
# initializes global variables
# Will be automatically called by init.game
repgames.startup = function() {
  if (exists("REP.CNT"))
    return();
  make.REP.CNT();
}

make.REP.CNT = function() {
  REP.CNT <<- list()
  
  # Tolerances
  REP.CNT$TOL <<- 10^(-9)
  REP.CNT$TOL.LY.get.upper.envelope <<- 10^(-8)
  REP.CNT$TOL.PM <<- 10^(-15)
  
    
  REP.CNT$PLOT.DURING.SOLVE <<- TRUE

  # From sk.glpk
  REP.CNT$LP.RETRY.COUNT <<- 0
  REP.CNT$LP.MAX.ITERATION <<- 10000
  REP.CNT$SOLVE.LP.WARNINGS <<- NULL
  REP.CNT$SINK_FILE <<- NULL
	REP.CNT$LP_PROBS <<- NULL
}


              