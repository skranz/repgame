# Interactively Solving Repeated Games: A Toolbox

# Author: Sebastian Kranz, Ulm University

This is an R package for numerically solving discounted infinitely repeated games in which players are risk neutral and can conduct voluntary monetary transfers in each round. It implements the algorithms developed in the article "Infinitely repeated games with public monitoring and monetary transfers" (JET, 2012) by Susanne Goldlücke and Sebastian Kranz.

## 1. Installation

### 1.1 Installing R, RTools and RStudio

First you need to install R, which is a very popular and powerful open source statistical programming language. You can download R for Windows, Max or Linux here:
  
  http://cran.r-project.org/

Note: If you have already installed R, you may want to update to the newest version by installing it again. 

If you use Windows, you also have to install RTools, which allow to compile C++ code, that is required for installing the repgame package from Github:

  https://cran.r-project.org/bin/windows/Rtools/

I recommend to additionally install RStudio, which is a great open source IDE for R:

http://rstudio.org/

### 1.2 Installing necessary R packages

You need to install several R packages from the internet. To do so, simply run in the R console the following code (you can use copy & paste):

```r

# Install CRAN Packages
pkgs = c("devtools","glpkAPI","slam","Rcpp","lattice")
for (pkg in pkgs) {
  if (!require(pkg, character.only=TRUE))
    install.packages(pkg)
}

# Install packages from GITHUB
library(devtools)
install_github(repo="skranz/restorepoint")
install_github(repo="skranz/rowmins")
install_github(repo="skranz/repgame")
```

# 2. Using the package

Here is a Tutorial from 2012 that explains an older version of the package in detail:

[Interactively Solving Repeated Games](http://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.160/pdf_dokumente/mitarbeiter/kranz/Interactively_Solving_Repeated_Games.pdf)

Installation instructions are outdated (use the ones above) and there where small changes in the package name, but most stuff hopefully still works. It is unfortunately not trivial to update the tutorial since my Lyx-R link is somehow broken. 

The subfolder examples contains updated versions of the tutorial's examples. Here is just the first example for an illustration:

```r
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
```

# Bugs and Issues

If you find a bug or have a suggestion, please send me an email (sebastian.kranz@uni-ulm.de) or file an issue for this github project.