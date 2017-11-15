plot.genarch <- function(x, min.vi = .5, main = NULL, cex.axis = 1, 
                         cex.names = 1, cex.main = 1, maxval = NULL, 
                         minval = NULL, col.ramp = NULL, ...){
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