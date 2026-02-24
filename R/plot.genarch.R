#' Plot Model-Averaged Line Cross Analysis Results
#'
#' S3 plot method for `"genarch"` objects produced by [LCA()]. Creates a
#' barplot of model-weighted average parameter estimates with unconditional
#' standard error bars. Bars are colored by variable importance, with hotter
#' colors indicating higher importance. Only composite genetic effects (CGEs)
#' with variable importance at or above `min.vi` are included.
#'
#' @param x A `"genarch"` object returned by [LCA()].
#' @param min.vi Numeric between 0 and 1. Minimum variable importance
#'   threshold for including a CGE in the plot. Default is `0.5`.
#' @param main Character. Main title for the plot. Default is
#'   `"Model Weighted Averages and Unconditional SE"`.
#' @param cex.axis Numeric. Expansion factor for the Y-axis labels. Default
#'   is `1`.
#' @param cex.names Numeric. Expansion factor for the CGE names on the
#'   X-axis. Default is `1`.
#' @param cex.main Numeric. Expansion factor for the main title. Default
#'   is `1`.
#' @param maxval Numeric. Custom maximum value for the Y-axis. If `NULL`,
#'   determined automatically from the data. Default is `NULL`.
#' @param minval Numeric. Custom minimum value for the Y-axis. If `NULL`,
#'   determined automatically from the data. Default is `NULL`.
#' @param col.ramp A vector of colors to use for the variable importance
#'   color scale. Default is `heat.colors(100)[100:1]`.
#' @param ... Additional arguments passed to [barplot()].
#'
#' @return Called for its side effect of producing a plot. Returns `NULL`
#'   invisibly.
#'
#' @seealso [LCA()] for performing the analysis, [VisModelSpace()] for
#'   model space visualization, [EvaluateModel()] for single model plots.
#'
#' @export
#' @examples
#' data(BanInf)
#' result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
#' plot(result)
#' plot(result, min.vi = 0.3)
plot.genarch <- function(x,
                         min.vi = .5,
                         main = NULL,
                         cex.axis = 1,
                         cex.names = 1,
                         cex.main = 1,
                         maxval = NULL,
                         minval = NULL,
                         col.ramp = NULL, ...){
  data <- x
  # if the user is getting rid of too much data lets stop and warn them
  if(min.vi > .5) cat("\nYou may want to reduce the min.vi argument to see other
                      important CGEs\n")
  #stop("select a lower min.vi to include more CGEs", call.=F)
  # reformat and keep only the data to be plotted
  x <- c(as.numeric(t(data$estimates)[-1, ]), as.numeric(data$varimp[,2]))
  x <- as.data.frame(matrix(x, length(x) / 3, 3))
  colnames(x) <- c("estimates", "SE", "vi")
  row.names(x) <- colnames(data$estimates)[-1]

  if(sum(x$vi >= min.vi) == 0){
    cat(paste("\nThe min.vi of", min.vi,
          "is larger than those observed",
          "in your dataset. Decrease the supplied",
          "value of min.vi to plot your results.\n"))
    stop()
  }
  data <- x[x$vi >= min.vi, ]
  # get some extra room
  par(mar=c(2, 2, 2, 6))
  # make colors for barplot
  if(is.null(col.ramp)) col.ramp <- heat.colors(100)[100:1]
  # set up graph max and min y axis if user does not provide
  if(is.null(maxval)){
    maxval <- max(as.numeric(data$estimates) + as.numeric(data$SE))
    if(maxval < 0) maxval <- 0
  }
  if(is.null(minval)){
    minval <- min(as.numeric(data$estimates) - as.numeric(data$SE))
    if(minval > 0) minval <- 0
  }
  # give the graph a title if the user fails to
  if(is.null(main)) main <- "Model Weighted Averages and Unconditional SE"
  mp <- barplot(as.numeric(data$estimates),
                names.arg=row.names(data),
                col=col.ramp[round(100 * as.numeric(data$vi))],
                ylim = c(minval - .4 * abs(minval), maxval + .4 * abs(maxval)),
                main = main,
                cex.axis=cex.axis,
                cex.names=cex.names,
                cex.main=cex.main)
  segments(mp, as.numeric(data$estimates) - as.numeric(data$SE),
           mp, as.numeric(data$estimates) + as.numeric(data$SE), lwd=2)
  # Now plot the horizontal bounds for the error bars
  # 1. The lower bar
  segments(mp - 0.1, as.numeric(data$estimates) - as.numeric(data$SE),
           mp + 0.1, as.numeric(data$estimates) - as.numeric(data$SE),
           lwd = 2)
  # 2. The upper bar
  segments(mp - 0.1, as.numeric(data$estimates) + as.numeric(data$SE),
           mp + 0.1, as.numeric(data$estimates) + as.numeric(data$SE),
           lwd = 2)
  # add a legend for the variable importance
  locs <- par("usr")
  color.legend(locs[2] ,  #xl
               locs[4]- ((locs[4]-locs[3]) * 0.5),  #yb
               locs[2] + (locs[2]*.05),    #xr
               locs[4],        #yt
               legend = c("0.00", "0.25", "0.50", "0.75","1"),
               rect.col = col.ramp,
               cex = 0.6,
               align="rb",
               gradient = "y")
  par(xpd = TRUE)
  text(x = (((locs[2]) + (locs[2] + (locs[2]*.05)))/2), y = locs[4],
       labels = "Variable\nImportance", cex = 0.5, pos = 3)
}
