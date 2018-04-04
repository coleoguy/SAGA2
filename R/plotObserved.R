plotObserved <- function(data, 
                         SE=T, 
                         pch=NULL, 
                         col=NULL, 
                         xlab=NULL, 
                         ylab=NULL, 
                         main=NULL,
                         env=FALSE){
  .pardefault <- par(no.readonly = T)
  
  # fill in some defaults
  if(is.null(pch)) pch <- 16
  if(is.null(col)) col <- "black"
  if(is.null(xlab)) xlab <- "% P1 Genome"
  if(is.null(ylab)) ylab <- "Phenotype Measure"
  
  # get the x and y
  x <- PrepareCmatrix(user.data = data, SCS = "XY", 
                      parental = "calc", drop.pars = NULL, 
                      messages=F)$p1a * 100
  y <- data$mean
  se <- data$SE
  
  # jitter the x values since they are often the same
  if(length(unique(x)) != length(x)){
     x <- jitter(x)
  }else{
     x <- x
  }
  
  if(SE==T){

    high <- max(y + se)
    low <- min(y - se)
    
    plot(x=x, y=y, 
         ylab=ylab, xlab=xlab, 
         xaxt="n", pch=pch, main=main, ylim=c(low,high), xlim=c(0,100))
    for(i in 1:length(x)){
      lines(x=rep(x[i], 2), y= c((y[i] + se[i]), (y[i] - se[i])))
    }
  }else{
    plot(x=x, y=y, ylab=ylab, xlab=xlab, xaxt="n", pch=pch, main=main)
  }
  axis(side=1,labels=c(0,50,100), at=c(0,50,100))
  abline(glm(y~x, weights = se), lty="dashed", col="blue")
  par(.pardefault)
}