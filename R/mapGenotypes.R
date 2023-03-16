##' Maximum posterior genotypes of unknown contributors
##'
##' @description For each marker, a ranked list of configurations of genotypes for some
##' or all unknown contributors is returned. The list contains all
##' configurations with posterior probability higher than some
##' specified \code{pmin}.
##'
##'\emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which
##' is intended as a service to enable users to try \pkg{DNAmixtures}
##' without purchasing a commercial licence for Hugin. When at all
##' possible, we strongly recommend the use of \pkg{DNAmixtures}
##' rather than this lite-version. See
##' \url{https://dnamixtures.r-forge.r-project.org/} for details on
##' both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of
##' \pkg{DNAmixtures}, note that computations are much less efficient
##' and that there are some differences in available functionality. Be
##' aware that the present documentation is copied from
##' \pkg{DNAmixtures} and thus may not accurately describe the
##' implementation of this lite-version.}
##'
##' @details Note that an error occurs if there are no configurations
##' with probability higher than \code{pmin}. In this case, try a
##' smaller \code{pmin}.
##'
##' The function makes use of
##' \code{\link[gRaven]{map.configurations}}, which localises the
##' configurations of highest high posterior probability by simulating
##' from the Bayesian networks until enough (at least mass
##' 1-\code{pmin}) of the state space has been explored -- the
##' computation time is thus dependent on how flat the posterior is,
##' and how small \code{pmin} is. The simulation is used only to locate
##' the relevant configurations; the computed probabilities are exact.
##'
##' @param mixture a DNA mixture
##' @param pmin A list of the minimum probability to consider for each marker.
##' @param U Optionally the indices of the unknown contributors of interest, specified as an integer vector.
##' @param markers Optionally, a subset of markers.
##' @param type It may be of interest to consider only the prediction of alleles in some subset of alleles. We allow
##' \describe{
##' \item{\code{"seen"}}{Consider only alleles that are seen in at least one EPG}
##' \item{\code{"all"}}{Consider the entire allelic range}
##' \item{\code{"unseen"}}{Consider only alleles that are not seen in any EPG (possibly redundant)}
##' }
##' @return A list, which for each marker contains the maximum
##' posterior configurations of allele counts (genotypes) above the specified probabilities \code{pmin}.
##' @export
##' @seealso \code{\link{summary.map.genotypes}}
##' @examples
##' data(MC18, USCaucasian, SGMplusDyes)
##' mix <- DNAmixture(list(MC18), k = 3, K = "K1", C = list(50), database = USCaucasian)
##' p <- mixpar(rho = list(30), eta = list(30), xi = list(0.07),
##'             phi = list(c(K1 = 0.7, U1 =0.2, U2 = 0.1)))
##' ## Inlude the peak height information
##' setPeakInfo(mix, p)
##' ## Marginally best genotypes for contributor U1
##' mpU1 <- map.genotypes(mix, pmin = 0.01, U = 1, type = "seen", markers = "D16S539")
##' summary(mpU1)
##' \donttest{
##' ## Jointly best genotypes for all unknown contributors
##' mp <- map.genotypes(mix, pmin = 0.01, type = "seen")
##' summary(mp) ## Profiles as genotypes rather than allelecounts
##' }
map.genotypes <- function(mixture, pmin,
                          U = seq_along(mixture$U), ## numeric index, not the names U1 etc
                          markers = mixture$markers, ## vector of marker-names
                          type = c("seen", "all", "unseen")
                          ){

  if (mixture$n.unknown == 0) stop("There are no unknown genotypes to predict")

  pmin <- rep(pmin, length.out = length(markers))
  if (!all(U < mixture$n.unknown)){"Invalid contributors"}
  type <- match.arg(type)

  getProbs <- function(m){
    ## extract peak heights
    H <- mixture$data[[m]][, seq_len(mixture$ntraces) + 1, drop = FALSE]

    ## Choose alleles
    ind <- switch(type,
                  seen = which(rowSums(H)>0),
                  unseen = which(rowSums(H)==0),
                  all = seq_len(nrow(H)))

    d <- map.configurations(mixture$domains[[m]],
                            attr(mixture$domains[[m]], "n")[ind, U], pmin = pmin[which(markers==m)])

    attr(d, "alleles") <- mixture$data[[m]]$allele[ind]
    d
  }
  out <- lapply(markers, getProbs)
  names(out) <- markers
  structure(out,
            U = U,
            pmin = pmin,
            ptotal = sapply(out, function(x)sum(x$Prob)),
            class = c("map.genotypes", "list")
            )
}

##' Summary of best genotypes
##'
##' The maximum posterior configurations of genotypes as returned by
##' \code{\link{map.genotypes}}, but in the format of pairs of alleles
##' rather than raw allelecounts. When a subset of alleles are
##' considered, then the value \code{NA} is used to denote an allele outside
##' this subset; in particular, if \code{type="seen"} is used in
##' \code{\link{map.genotypes}}, then \code{NA} corresponds to the
##' allele having dropped out in all of the mixtures included in the
##' model.
##'
##' @method summary map.genotypes
##' @param object An object returned by \code{\link{map.genotypes}}.
##' @param ... arguments passed on to other methods.
##' @return A \code{data.frame} with two columns per unknown
##' contributor under consideration and a column \code{Prob}
##' containing the probability of the configuration.
##' @seealso \code{\link{map.genotypes}}
##' @author Therese Graversen
##' @export
summary.map.genotypes <- function(object, ...){

  extractGenotype <- function(dat){
    alleles <- attr(dat, "alleles")
    vec <- rep(seq_along(attr(object, "U")), each = length(alleles))
    mat <- matrix(NA, nrow = nrow(dat), ncol = 2*length(attr(object, "U")))
    for (i in 1:nrow(dat)){
      r <- split(head(as.numeric(dat[i,]), -1), vec)
      for (j in seq_along(attr(object, "U"))){
        mat[i, 2*(j-1) + 1:2] <- rep(c(alleles, NA), times = c(r[[j]], 2-sum(r[[j]])))
      }
    }
    mat <- as.data.frame(mat)
    names(mat) <-  paste0("U", rep(attr(object, "U"), each = 2), ".",
                          rep(1:2, times = length(attr(object, "U"))))
    mat$Prob <- dat$Prob
    mat
  }
  out <- lapply(object, extractGenotype)

  structure(out,
            class = c("summary.map.genotypes", "list"),
            U = attr(object, "U"),
            pmin = attr(object, "pmin"),
            ptotal = attr(object, "ptotal")
            )
}

##' @rdname summary.map.genotypes
##' @method print summary.map.genotypes
##' @param x An object of class \code{"summary.map.genotypes"}, typically returned by \code{summary.map.genotypes}.
##' @param markers Optionally, a subset of markers to print
##' @export
print.summary.map.genotypes <- function(x, markers = names(x), ...){
  for (i in markers){
    cat("\n", i, ":\n", sep = "")
    print(as.data.frame(do.call(cbind, c(lapply(x[[i]][,!names(x[[i]])=="Prob", drop=FALSE], format, drop0trailing = TRUE),
                                         list(Prob = format(x[[i]]$Prob, digits=4))))), right = FALSE, quote=FALSE,
          print.gap = 3L, ...)
    cat("\nTotal probability:", format(attr(x, "ptotal")[i], digits = 4), "\n\n")
  }
}


