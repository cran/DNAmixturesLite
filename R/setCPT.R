##' Set conditional probability tables for all auxiliary variables
##'
##' @description A wrapper-function for the case where the conditional probability
##' tables are to be set in a DNA mixture model for all the auxiliary variables
##' (O, D and Q, at all markers and all alleles).
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @details
##' The conditional probability tables are set according to a list of
##' specified parameters and the observed peak heights. The function
##' returns evidence for future use on the nodes O.
##'
##' @export
##' @param mixture A DNA mixture
##' @param pars A list of parameters
##' @param markers Optionally, a list of markers for which to set the tables
##' @return List containing evidence for conditioning on peak heights.
##' @seealso For further details, see in particular \code{\link{setCPT.O}}.
setCPT <- function(mixture, pars, markers = mixture$markers){

  ## extract relevant stuff from mixture and pass on
  C <- mixture$C
  k <- mixture$k
  K <- mixture$K
  U <- mixture$U
  n.unknown <- mixture$n.unknown
  observed <- mixture$observed

  data <- mixture$data
  domains <- mixture$domains

  one.marker <- function(m){
    domain <- domains[[m]]
    d <- data[[m]]
    n_K <- subset(d, select = K)
    gs <- d$gets_stutter
    cs <- d$can_stutter
    sf <- d$stutter.from

    ev <- list()
    for (r in observed[[m]]){
      rho <- pars[[r,"rho"]]
      xi <- pars[[r,"xi"]]
      eta <- pars[[r,"eta"]]
      phi <- pars[[r,"phi"]][c(U, K)]

      setCPT.D(domain, rho, xi, eta, phi, C[[r]], n.unknown,
                n_K, attr(domain, "D")[[r]], gs, cs, sf)

      setCPT.Q(domain, rho, xi, eta, phi, d[,r+1], n.unknown,
                n_K, attr(domain, "Q")[[r]], gs, cs, sf)

      ev[[r]] <- setCPT.O(domain, rho, xi, eta, phi, d[,r+1], C[[r]], n.unknown,
                           n_K, attr(domain, "O")[[r]], gs, cs, sf)

    }
    ev
  }
  evidence <- lapply(markers, one.marker)
  names(evidence) <- markers
  evidence
}

##' Include or exclude peak information in the model
##'
##' @description The function \code{setPeakInfo} is used for including in the model
##' either the full peak height information or only information about
##' presence and absence of peaks. After a call to \code{setPeakInfo},
##' the Bayesian networks in \command{mixture$domains} will represent
##' the conditional distribution given the specified peak
##' information. For the reverse functionality, i.e. exclusion of any
##' such peak information, use \code{removePeakInfo}.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @details The function \code{setPeakInfo} sets conditional
##' probability tables using the specified model parameters, and
##' propagates suitable evidence using either observed peak heights or
##' the discrete presence/absence observations of peaks. Any
##' previously entered or propagated evidence on nodes \code{O} and
##' \code{D} will be retracted in this process.
##'
##' The function \code{removePeakInfo} retracts all evidence on nodes
##' \code{O} and \code{D}; the conditional probability tables are left
##' unchanged.
##'
##' @param mixture A \code{\link{DNAmixture}}
##' @param pars A \code{\link{mixpar}} model parameter
##' @param presence.only Default is FALSE, which means that the full
##' peak height information is taken into consideration. Set to TRUE,
##' which will include only information on the presence and absence of peaks.
##' @return invisibly \code{NULL}
##' @seealso \code{\link{setCPT}}. For use of the Bayesian networks, see \code{\link{map.genotypes}}.
##' @export
setPeakInfo <- function(mixture, pars, presence.only = FALSE){
  ## Set conditional probability tables for O, D, Q.
  evidence <- setCPT(mixture, pars)

  for (m in mixture$markers){
    D <- attr(mixture$domains[[m]], "D")
    O <- attr(mixture$domains[[m]], "O")
    retract(mixture$domains[[m]], unlist(c(D, O)))

    for (r in mixture$observed[[m]]){
      for (a in seq_along(O[[r]])){
        if (presence.only){
          set.finding(mixture$domains[[m]], D[[r]][a], mixture$data[[m]][a, r+1] <= 0)
        }
        else {
          set.finding(mixture$domains[[m]], O[[r]][a], evidence[[m]][[r]][,a])
        }
      }
    }
    propagate(mixture$domains[[m]])
  }
  invisible(NULL)
}

##' @export
##' @rdname setPeakInfo
removePeakInfo <- function(mixture){

  for (m in mixture$markers){
    D <- attr(mixture$domains[[m]], "D")
    O <- attr(mixture$domains[[m]], "O")
    retract(mixture$domains[[m]], unlist(c(D, O)))
    propagate(mixture$domains[[m]])
  }
  invisible(NULL)
}
