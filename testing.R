# figure out if how we can add vectors of parental strains
# this would allow us to handle pooled strains like in the 
# mojavensis line


# Load our new function
source("R/AnalyzeModels.R")
source("R/support.R")
source("R/PrepareCmatrix.R")
source("R/LCA.R")
source("R/plot.genarch.R")
library(plotrix) #color.legend
library(gtools)  #combinations

# This is the example of user data
data1 <- read.csv("emp.data/ban.osa.29.csv", as.is = T)
data2 <- read.csv("emp.data/ban.osa.35.csv", as.is = T)
data3 <- read.csv("emp.data/ban.osa.comb.csv", as.is=T)
res1 <- LCA(data = data1, SCS = "XY", parental = "calc", model.sum = .95, env = T, drop.pars=c("Med"), max.pars = 5)
res1 <- LCA(data = data1, SCS = "XY", parental = "calc", model.sum = .95, max.pars = 7)
res2 <- LCA(data = data2, SCS = "XY", parental = "calc", model.sum = .95, max.pars = 7)
res3 <- LCA(data = data3, SCS = "XY", parental = "calc", model.sum = .95, max.pars = 5)
par(mfcol=c(3,1))
plot(res1,min.vi = .5)
plot(res2,min.vi = .5)
plot(res3,min.vi = .5)

# This secrion is for working on PrepareCmatrix
user.data <- data3
SCS<-"SCS"
parental<-"calc"
drop.pars<- "Meo"
env<-T

# arguments for AnalyzeModels
data=
SCS=
Cmatrix=
crosses=
parental=
env=
model.sum=
max.models=
max.pars=
  
# arguments for PrepareCmatrix
user.data=
SCS=
parental=
drop.pars=
env=
  
# arguments for LCA
SCS = "XY"
parental = "calc"
env = FALSE
model.sum = .95
max.models = 300000
drop.pars = NULL
max.pars = NULL
Cmatrix = NULL
      
      
  VisModelSpace(res3)
  
