#' Evaluate a Single Model from Line Cross Analysis
#'
#' Plots the parameter estimates and conditional standard errors for a single
#' model selected from the results of [LCA()]. Also returns the estimates,
#' a chi-squared goodness-of-fit p-value, and a Bonferroni-corrected alpha
#' level.
#'
#' @param data A `"genarch"` object returned by [LCA()] with `ret.all = TRUE`.
#' @param model Integer. The index of the model to evaluate, corresponding to
#'   the position in the `all.models` element of the `"genarch"` object.
#' @param cex.axis Numeric. Expansion factor for the Y-axis labels. Default
#'   is `1`.
#' @param cex.names Numeric. Expansion factor for the CGE names on the
#'   X-axis. Default is `1`.
#' @param cex.main Numeric. Expansion factor for the main title. Default
#'   is `1`.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{Estiamtes}{A 2-row matrix with parameter estimates (row 1) and
#'     conditional standard errors (row 2) for each CGE in the model.}
#'   \item{P.val}{Numeric. P-value from a chi-squared goodness-of-fit test
#'     of the model's residual deviance.}
#'   \item{Bonf.Corr.Alpha}{Numeric. Bonferroni-corrected significance level
#'     based on the total number of models evaluated.}
#' }
#'
#' @seealso [LCA()] for performing the analysis, [VisModelSpace()] for
#'   identifying models of interest.
#'
#' @export
#' @examples
#' data(BanInf)
#' result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
#' # Evaluate the first model
#' EvaluateModel(result, model = 1)
EvaluateModel <-
function(data, model, cex.axis=1, cex.names=1, cex.main=1){
  .pardefault <- par(no.readonly = T)
  result <- data[[1]][[model]]
  pars <- c("Mean", strsplit(names(data[[1]])[model], split = ", ")[[1]])
  par.est <- summary(result)$coefficients[, 1]
  se.est <- summary(result)$coefficients[, 2]
  max.val <- max(par.est + se.est)
  min.val <- min(par.est - se.est)
  if(max.val < 0){
    max.val <- 0
  }
  if(min.val > 0){
    min.val <- 0
  }
  mp <- barplot(par.est, names.arg = pars,
                main = "Single Model Means and Cond. SE",
                ylim = c(min.val - 0.2, max.val + 0.2), cex.axis=cex.axis,
                cex.names=cex.names, cex.main=cex.main)

  high.se <- par.est + se.est
  low.se <- par.est - se.est
  segments(mp, high.se, mp, low.se, lwd = 3)
  segments(mp - 0.1, high.se, mp + 0.1, high.se, lwd = 3)
  segments(mp - 0.1, low.se, mp + 0.1, low.se, lwd = 3)
  results <- matrix(, 2, length(pars))
  results[1, ] <- par.est
  results[2, ] <- se.est
  colnames(results) <- pars
  row.names(results) <- c("Estimate", "Cond. SE")
  model.p <- 1 - pchisq(result$deviance, result$df.residual, lower.tail = T)
  mult.test <- 0.05 / length(data)
  final.results <- list()
  final.results[[1]] <- results
  final.results[[2]] <- model.p
  final.results[[3]] <- mult.test
  names(final.results) <- c("Estiamtes", "P.val", "Bonf.Corr.Alpha")
  return(final.results)
  par(.pardefault)
}
