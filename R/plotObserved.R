#' Plot Observed Line Cross Data
#'
#' Creates a traditional line cross analysis plot showing observed cohort
#' phenotype means as a function of percentage of parent 1 genome. Overlays
#' the expectation from a simple additive model as a dashed blue line.
#' Optionally includes standard error bars.
#'
#' @param data A data frame in the format required by [LCA()], with columns
#'   `cross`, `mean`, `SE`, `sex`, `enviro`, `sire`, and `dam`.
#' @param SE Logical. Whether to include standard error bars. Default is
#'   `TRUE`.
#' @param pch Integer or character. Plotting symbol for observed means. If
#'   shorter than the number of cohorts, it will be recycled. Default is `16`.
#' @param col Character. Color for plotting symbols. If shorter than the
#'   number of cohorts, it will be recycled. Default is `"black"`.
#' @param xlab Character. Label for the X-axis. Default is
#'   `"% P1 Genome"`.
#' @param ylab Character. Label for the Y-axis. Default is
#'   `"Phenotype Measure"`.
#' @param main Character. Main title for the plot. Default is `NULL` (no
#'   title).
#' @param env Logical. Whether an environmental effect is present in the
#'   data. Default is `FALSE`.
#'
#' @return Called for its side effect of producing a plot. Returns `NULL`
#'   invisibly.
#'
#' @seealso [LCA()] for performing the analysis, [plot.genarch()] for
#'   plotting model-averaged results.
#'
#' @export
#' @examples
#' data(BanInf)
#' plotObserved(BanInf)
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
                      keep.pars = NULL, messages=F, env=F)$Aa
  x <- (x-min(x)) / (max(x)-min(x)) * 100
  y <- data$mean
  se <- data$SE

  # jitter the x values since they are often the same
  if(length(unique(x)) != length(x)){
     x <- x + rnorm(mean=0, sd=.8, n=length(x))
  }else{
     x <- x
  }

  if(SE==T){

    high <- max(y + se)
    low <- min(y - se)

    plot(x=x, y=y, col=rgb(.1,.1,.1,.5),
         ylab=ylab, xlab=xlab,
         xaxt="n", pch=pch, main=main, ylim=c(low,high), xlim=c(0,100))
    for(i in 1:length(x)){
      lines(x=rep(x[i], 2), y= c((y[i] + se[i]), (y[i] - se[i])), col=rgb(.1,.1,.1,.5))
    }
  }else{
    plot(x=x, y=y, ylab=ylab, xlab=xlab, col=rgb(.1,.1,.1,.5), xaxt="n", pch=pch, main=main)
  }
  axis(side=1,labels=c(0,50,100), at=c(0,50,100))
  abline(glm(y~x, weights = se), lty="dashed", col="blue")
  par(.pardefault)
}
