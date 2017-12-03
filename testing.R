library(SAGA2)

dat1 <- read.csv("inst/ba.sq.26.csv", as.is=T)
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











