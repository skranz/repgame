# Interactively Solving Repeated Games: A Toolbox

# Author: Sebastian Kranz, Ulm University

## 1. Installing neccessary software

### 1.1 Installing R and RStudio

First you need to install R, which is a very popular and powerful open source statistical programming language. You can download R for Windows, Max or Linux here:
  
  http://cran.r-project.org/

Note: If you have already installed R, you may want to update to the newest version by installing it again. 

I recommend to additionally install RStudio, which is a great open source IDE for R:

http://rstudio.org/

### 1.2 Installing necessary R packages

You need to install several R packages from the internet. To do so, simply run in the R console the following code (you can use copy & paste):
```{r eval=FALSE}

# Install CRAN Packages
pkgs = c("devtools","glpkAPI","slam","RCPP","lattice")
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
