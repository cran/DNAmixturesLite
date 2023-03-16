##' Calculate prequential scores
##'
##' @description \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @param mixture A \code{DNAmixture} model object.
##' @param pars A \code{mixpar} model parameter.
##' @param markers An ordering of the markers, possibly a subset of the markers only.
##' @param by.allele Should conditioning be done allele-wise (\code{TRUE}) or EPG-wise (\code{FALSE}).
##' Defaults to \code{TRUE}. For details, see \code{\link{predict}}.
##' @return A data.frame, which contains the output from
##' \code{\link{predict}} as well as columns \code{Y}, \code{EY}, \code{VY}, corresponding
##' to the log-score and its mean and variance. Finally the variable
##' \code{score} is added, which is the normalised cumulative log-score.
##' @author Therese Graversen
##' @export
##' @examples
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50,50), k = 3, K = c("K1", "K3"),
##'                   database = USCaucasian)
##' p <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08,0.08),
##'             phi = list(c(U1 = 0.1, K3 = 0.2, K1 = 0.7), c(U1 = 0.1, K3 = 0.2, K1 = 0.7)))
##' preq <- prequential.score(mix, pars = p)
##' plot(preq, col = preq$trace)
##' ## Annotate using repeat numbers
##' text(preq$score, labels = preq$allele, pos = c(1,3), col = preq$trace, cex = 0.6)
prequential.score <- function(mixture, pars, markers = mixture$markers, by.allele = TRUE){

  ## We need probabilities of not seeing a peak conditionally on previous alleles
  ## p_a = P(Z_a > 0| Z_b <= z_b for b < a)
  stopifnot(class(mixture) == "DNAmixture")
  out <- predict(mixture, pars, dist = "prequential", by.allele = by.allele, markers = markers)

  ## Defining log0 = 0, will only be used as 0*log0
  logzero <- function(x)ifelse(x == 0, 0, log(x))

  one.marker <- function(d){
    d$Y <- - log(ifelse(d$height > 0, d$seen, d$unseen))
    d$EY <- - d$seen*logzero(d$seen) - (d$unseen)*logzero(d$unseen)
    d$VY <- d$seen*logzero(d$seen)^2 + d$unseen*logzero(d$unseen)^2 - d$EY^2
    d
  }
  out <- lapply(out, one.marker)
  out <- do.call(rbind, out)
  rownames(out) <- NULL

  out$score <- ifelse(cumsum(out$VY) == 0, 0, cumsum(out$Y - out$EY))
  class(out) <- c("prequential.score", "data.frame")
  out
}


##' @rdname prequential.score
##' @method plot prequential.score
##' @export
##' @param x A data.frame containing at least variables \code{marker} and \code{score}.
##' @param normalise Should the prequential score be normalised? Defaults to \code{FALSE}.
##' @param ylab Label for the y-axis.
##' @param ylim Range for the y-axis.
##' @param ... Additional arguments to be passed on to \code{plot}.
plot.prequential.score <- function(x, normalise = FALSE, ylab = NULL, ylim = NULL, ...){

  markers <- unique(x$marker)
  M <- match(markers, x$marker)

  if (normalise){
    if (is.null(ylab)) ylab <- "Normalised cumulative score"
    plot(x$score,
         type = "n", xaxt = "n",
         xlab = "", ylab = ylab, ...)
    abline(v = M, lty = "22", col = "grey")
    abline(h = qnorm(c(0.95,0.99)))
    points(x$score/sqrt(cumsum(x$VY)), type = "b", ...)
    axis(1, at = M + c(diff(M)/2, (nrow(x) - M[length(M)])/2), labels = markers,
         tick = FALSE, ...)
    return(invisible(x))
  }
  else{
    if (is.null(ylab)) ylab <- "Cumulative score"
    if (is.null(ylim)) ylim <- range(x$score, qnorm(0.99)*sqrt(sum(x$VY)))
    plot(x$score,
         type = "n", xaxt = "n",
         xlab = "", ylab = ylab, ylim = ylim, ...)
    abline(v = M, lty = "22", col = "grey")
    lines(qnorm(0.95)*sqrt(cumsum(x$VY)), lty = "43")
    lines(qnorm(0.99)*sqrt(cumsum(x$VY)), lty = "22")
    points(x$score, type = "b", ...)
    axis(1, at = M + c(diff(M)/2, (nrow(x) - M[length(M)])/2), labels = markers,
         tick = FALSE, ...)
    return(invisible(x))
  }

}
