library(SAGA2)

dat1 <- read.csv("ba.sq.26.csv", as.is=T)
dat2 <- read.csv("inst/ba.sq.29.csv", as.is=T)
dat3 <- read.csv("inst/ba.sq.35.csv", as.is=T)
dat4 <- read.csv("inst/ba.sq.comb.csv", as.is=T)

dat5 <- read.csv("inst/bho.sq.26.csv", as.is=T)
dat6 <- read.csv("inst/bho.sq.29.csv", as.is=T)
dat7 <- read.csv("inst/bho.sq.35.csv", as.is=T)
dat8 <- read.csv("inst/bho.sq.comb.csv", as.is=T)

dat9  <- read.csv("inst/da.sq.26.csv", as.is=T)
dat10 <- read.csv("inst/da.sq.29.csv", as.is=T)
dat11 <- read.csv("inst/da.sq.35.csv", as.is=T)
dat12 <- read.csv("inst/da.sq.comb.csv", as.is=T)

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
plot(res1, main=paste("ba x sq 26C:", length(res1$best.models), "models"),cex.names=.6)
plot(res2, main=paste("ba x sq 29C", length(res2$best.models), "models"),cex.names=.6)
plot(res3, main=paste("ba x sq 35C", length(res3$best.models), "models"),cex.names=.6)
plot(res4, main=paste("ba x sq all data", length(res4$best.models), "models"),cex.names=.6)
plot(res5, main=paste("bho x sq 26C", length(res5$best.models), "models"),cex.names=.6)
plot(res6, main=paste("bho x sq 29C", length(res6$best.models), "models"),cex.names=.6)
plot(res7, main=paste("bho x sq 35C", length(res7$best.models), "models"),cex.names=.6)
plot(res8, main=paste("bho x sq all data", length(res8$best.models), "models"),cex.names=.6)
plot(res9, main=paste("da x sq 26C", length(res9$best.models), "models"),cex.names=.6)
plot(res10, main=paste("da x sq 29C", length(res10$best.models), "models"),cex.names=.6)
plot(res11, main=paste("da x sq 35C", length(res11$best.models), "models"),cex.names=.6)
plot(res12, main=paste("da x sq all data", length(res12$best.models), "models"),cex.names=.6)




# old SAGA
Aa, Ad, Xa, Xd, Ca, Mea, Med, AaAa, AaAd, AdAd, XaAa, CaAa, XaAd, XdAa, XdAd, CaAd, CaXa, CaXd
Aa, Ad, Xa, Xd, Ca, Mea, Med, AaAa, AaAd, AdAd, AaXa, AaCa, AdXa, AaXd, AdXd, AdCa, XaCa, XdCa, XaXa, XaXd, XdXd 
# lets compare 2 versions w/o preallocation
dat1 <- read.csv("silene.ovules.csv", as.is=T)
s.time <- Sys.time()
res4 <- LCA(dat1, max.pars=3)
e.time <- Sys.time()
e.time-s.time
# 10.419
# 10.498
# 10.69
# 10.79
# 11.76
# 10.78
# 10.89
# 10.47


library(devtools)
install_github('coleoguy/SAGA2')
library(SAGA2)
dat1 <- read.csv("silene.ovules.csv", as.is=T)
s.time <- Sys.time()
res4 <- LCA(dat1, max.pars=4)
e.time <- Sys.time()
e.time-s.time
# 10.19
# 11.24
# 10.94
# 9.962
#