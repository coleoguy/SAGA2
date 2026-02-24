#' Visualize Model Space
#'
#' Plots a grid of colored squares representing every model evaluated during
#' line cross analysis. The color of each square reflects its Akaike weight,
#' with hotter colors indicating higher support. Models are arranged from
#' simplest (bottom-left) to most complex (top-right).
#'
#' @param data A `"genarch"` object returned by [LCA()].
#' @param cex.u Numeric. Expansion factor for the model squares. Default
#'   is `3`.
#' @param cex.mtext Numeric. Expansion factor for axis labels. Default
#'   is `1`.
#' @param cex.leg Numeric. Expansion factor for the legend text. Default
#'   is `0.8`.
#' @param cex.mod Numeric. Expansion factor for model number labels within
#'   each square. Default is `0.35`.
#'
#' @return Called for its side effect of producing a plot. Returns `NULL`
#'   invisibly.
#'
#' @seealso [LCA()] for performing the analysis, [EvaluateModel()] for
#'   investigating individual models.
#'
#' @export
#' @examples
#' data(BanInf)
#' result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
#' VisModelSpace(result)
VisModelSpace <- function(data,
                          cex.u = 3,
                          cex.mtext=1,
                          cex.leg=.8,
                          cex.mod=.35){
  .pardefault <- par(no.readonly = T)
  models <- data[[1]] #mod.results
  data <- data[[5]]   #daicc
  colors <- heat.colors(1101)
  colors <- colors[1001:1]
  waic <- (exp(-.5 * data) / (sum(exp(-.5 * data))))
  foobar <- waic/max(waic)
  mod.heat <- round(1000 * foobar) + 1
  foo <- round(sqrt(length(data)) + 1)
  yvals <- rep(1:foo, each=foo)[1:length(data)]
  xvals <- rep(1:foo, times=foo)[1:length(data)]
  par(mar=c(2, 2, 2, 7))

  ## add the increasing level of model complexity

  plot(x=xvals, y=yvals, col=colors[mod.heat], pch=15,
       xaxt="n", yaxt="n", cex = cex.u,
       )
  text(x=xvals, y=yvals, labels=1:length(xvals), cex=cex.mod)
  x="Increasing Model Complexity"
  y=""
  foobar <- substitute(x %->% y)
  mtext(text=foobar, line=0.45, side=1, cex=cex.mtext)
  mtext(text=foobar, line=0, side=2, cex=cex.mtext)
  abbi <- 1:length(data)
  locs <- par("usr")
  color.legend(locs[2] + (locs[2]*.01),  #xl
               locs[4]- ((locs[4]-locs[3]) * 0.5),  #yb
               locs[2] + (locs[2]*.05),    #xr
               locs[4],        #yt
               legend = c("0.00", round(max(waic), digits=3)),
               rect.col = colors[1:length(colors)],
               cex = cex.leg,
               align = "rb", gradient = "y")
  par(xpd = TRUE)
  text(x=(((locs[2] + 1.5) + (locs[2] + (locs[2]*.05)))/2),
       y=(locs[4]+.2), labels="Akaike Weights", pos=3, cex=cex.leg)
  par(.pardefault)
}
