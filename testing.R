library(SAGA2)
dat1 <- read.csv("inst/pm.sin.29.csv", as.is=T)
dat1$enviro <- NA
res1 <- LCA(dat1)
plot(res1)


dat1 <- read.csv("ba.sq.26.csv", as.is=T)
dat2 <- read.csv("ba.sq.29.csv", as.is=T)
dat3 <- read.csv("ba.sq.35.csv", as.is=T)
dat4 <- read.csv("ba.sq.comb.csv", as.is=T)

dat5 <- read.csv("bho.sq.26.csv", as.is=T)
dat6 <- read.csv("bho.sq.29.csv", as.is=T)
dat7 <- read.csv("bho.sq.35.csv", as.is=T)
dat8 <- read.csv("bho.sq.comb.csv", as.is=T)

dat9  <- read.csv("da.sq.26.csv", as.is=T)
dat10 <- read.csv("da.sq.29.csv", as.is=T)
dat11 <- read.csv("da.sq.35.csv", as.is=T)
dat12 <- read.csv("da.sq.comb.csv", as.is=T)

res1 <- LCA(dat1)
res2 <- LCA(dat2)
res3 <- LCA(dat3)
res4 <- LCA(dat4)

res5 <- LCA(dat5)
res6 <- LCA(dat6)
res7 <- LCA(dat7)
res8 <- LCA(dat8)

res9  <- LCA(dat9)
res10 <- LCA(dat10)
res11 <- LCA(dat11)
res12 <- LCA(dat12)


par(mfcol=c(4,3))
plot(res1, main=paste("ba x sq 26C:", length(res1$best.models), "models"))
plot(res2, main=paste("ba x sq 29C", length(res2$best.models), "models"))
plot(res3, main=paste("ba x sq 35C", length(res3$best.models), "models"))
plot(res4, main=paste("ba x sq all data", length(res4$best.models), "models"))
plot(res5, main=paste("bho x sq 26C", length(res5$best.models), "models"))
plot(res6, main=paste("bho x sq 29C", length(res6$best.models), "models"))
plot(res7, main=paste("bho x sq 35C", length(res7$best.models), "models"))
plot(res8, main=paste("bho x sq all data", length(res8$best.models), "models"))
plot(res9, main=paste("da x sq 26C", length(res9$best.models), "models"))
plot(res10, main=paste("da x sq 29C", length(res10$best.models), "models"))
plot(res11, main=paste("da x sq 35C", length(res11$best.models), "models"))
plot(res12, main=paste("da x sq all data", length(res12$best.models), "models"))











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
  
