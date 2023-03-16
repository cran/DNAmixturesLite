##' Shape parameters for a mixture with known contributors only
##'
##' @description The function is mainly for internal use. It calculates the shape
##' parameters, or contribution to the shape parameter, that
##' correspond to the known contributors in the mixture.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##'
##' @param mixture A \code{DNAmixture}.
##' @param pars Model parameters.
##' @param allelecounts This argument is possibly redundant.
##' @return For each marker a matrix of shape-parameters.
##' @author Therese Graversen
##' @export
##' @seealso predict
##' @examples
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50,50), k = 3, K = c("K1", "K3", "K2"),
##' database = USCaucasian)
##' p <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08,0.08),
##'             phi = list(c(K2 = 0.1, K3 = 0.2, K1 = 0.7), c(K2 = 0.1, K3 = 0.2, K1 = 0.7)))
##' sh <- getShapes(mix, p)
##' sh$VWA
getShapes <- function(mixture, pars, allelecounts = NULL){

  dat <- mixture$data
  if (!missing(allelecounts)){
    dat <- lapply(dat, merge(data, allelecounts, all = TRUE))
  }

  rho <- pars[,"rho"]
  xi <- pars[,"xi"]
  ## order phi according to contributors (done trace-wise)
  phi <- lapply(pars[,"phi"], "[", mixture$K)

  one.marker <- function(m){
    d <- dat[[m]]
    ## The data is ordered according to phi.
    ns <- as.matrix(subset(d, select = mixture$K))
    shapes <- matrix(nrow = mixture$ntraces, ncol = nrow(d))

    ## non-observed combinations of marker and DNA sample has NA
    for (r in mixture$observed[[m]]){
      for (a in seq_len(nrow(d))){
        shapes[r,a] <- rho[[r]] * (1 - xi[[r]]*d$can_stutter[a]) * (ns[a,] %*% phi[[r]])
        if (d$gets_stutter[a]){
          st <- d$stutter.from[a]
          shapes[r,a] <- shapes[r,a] + rho[[r]] * xi[[r]] * (ns[st,] %*% phi[[r]])
        }
      }
    }
    shapes
  }
  ans <- lapply(names(dat), one.marker)
  names(ans) <- names(dat)
  ans
}

