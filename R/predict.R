##' Various probabilities in a fitted DNA mixture model
##'
##' @description
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @details For a mixture with unknown contributors, the
##' probabilities are computed with respect to one of three
##' distributions. Let \code{height} be the matrix of peak heights
##' with columns \code{height1, \ldots, heightR}. For a peak at allele \code{a} in the mixture \code{r}, the three choices of distributions are
##' \describe{
##' \item{\code{"joint"}}{Default. No conditioning on observed peak heights.}
##' \item{\code{"conditional"}}{Conditional on \code{height[-a, -r]}, i.e. on heights for all peaks, except the one under consideration.}
##' \item{\code{"prequential"}}{Conditional on \code{height[1:(a-1), 1:(r-1)]}, i.e. on heights for all peaks "before" the peak under consideration (see argument \code{by.allele} for details).}
##' }
##' If all contributors are known, the three distributions are the same
##' due to independence of the peak heights.
##'
##' @param object A \code{\link{DNAmixture}} object
##' @param pars Array of model parameters
##' @param dist One of "joint", "conditional", and "prequential". If
##' there are only known contributors, these are all the same since, under the model, peak
##' heights are condtionally independent given profiles of the
##' contributors.
##' @param markers The set of markers of interest
##' @param by.allele If \code{dist = "prequential"} then the order in
##' which we condition on mixtures and alleles matters. \code{by.allele
##' = TRUE} will proceed through alleles in increasing repeat number,
##' and for each allele condition on one mixture at the time. If \code{FALSE},
##' the conditioning is done by mixtures and then alleles within these.
##' @param initialize By default \code{predict} removes all entered
##' evidence from the networks in \code{object}. Setting
##' \code{initialize = FALSE} should be done with care, and it is up
##' to the user to ensure that the returned probabilities are meaningful.
##' @param ... Not used
##'
##' @return A list with one data.frame per marker containing various probabilities for
##' diagnostics
##' \item{unseen}{The probability of not seeing a peak, i.e. no peak or a peak falling below the threshold}
##' \item{seen}{The probability of seeing the allele}
##' \item{smaller}{The probability of seeing a smaller peak than the
##' one observed}
##' \item{larger}{The probability of seeing a larger
##' peak than the one observed}
##' @author Therese Graversen
##' @method predict DNAmixture
##' @export
##' @examples
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50,50), k = 3, K = c("K1", "K3", "K2"),
##' database = USCaucasian)
##' p <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08,0.08),
##'             phi = list(c(K2 = 0.1, K3 = 0.2, K1 = 0.7), c(K2 = 0.1, K3 = 0.2, K1 = 0.7)))
##' pred <- predict(mix, p)
##' pred$VWA
predict.DNAmixture <- function(object, pars, dist = c("joint", "conditional", "prequential"),
                               markers = object$markers, by.allele = TRUE,
                               initialize = TRUE, ...){

  observed <- object$observed

  if (object$n.unknown == 0){
    dist <- "known"
    all.shapes <- getShapes(object, pars)
    eta <- unlist(pars[,"eta"])
    C <- unlist(object$C)
  }
  else {
    dist <- match.arg(dist)
    ## Remove all propagated evidence in the network
    if (initialize) lapply(object$domains[markers], initialize.domain)
    ## Note that tables are set on the global mixture object (and all markers).
    evidence <- setCPT(object, pars, markers = markers)
    ## Now we need to propagate before querying for beliefs,
    ## this is done below.
  }

  f.known <- function(m){
    traces <- observed[[m]]
    ## Inner loop will be traces, and outer loop will be alleles
    heights <- as.vector(t(object$data[[m]][,traces+1])) ## a column per trace
    shapes <- as.vector(all.shapes[[m]][traces,])        ## a row per trace, with NA-s for non-observed
    alleles <- object$data[[m]][,1]

    cbind(
      rep(alleles, each = length(traces)),
      pgamma(C[traces], shape = shapes, scale = eta[traces], lower.tail = FALSE),
      pgamma(C[traces], shape = shapes, scale = eta[traces]),
      pgamma(heights, shape = shapes, scale = eta[traces], lower.tail = FALSE),
      pgamma(heights, shape = shapes, scale = eta[traces]),
      rep(traces, times = length(alleles)),
      heights
      )
  }

  f.joint <- function(m){
    domain <- object$domains[[m]]
    D <- attr(domain, "D")
    Q <- attr(domain, "Q")
    dat <- object$data[[m]]

    propagate(domain)

    out <- list()
    for (r in observed[[m]]){
      ds <- D[[r]]
      qs <- Q[[r]]
      out[[r]] <- do.call("rbind", lapply(seq_along(ds), function(a){
        c(dat$allele[a], get.belief(domain, ds[a]), get.belief(domain, qs[a]), r, dat[a,r+1])
      }))
    }
    do.call("rbind", out)
  }

  f.conditional <- function(m){
    domain <- object$domains[[m]]
    D <- attr(domain, "D")
    Q <- attr(domain, "Q")
    O <- attr(domain, "O")
    dat <- object$data[[m]]

    ## Enter evidence on all peak heights
    for (r in observed[[m]]){
      o <- O[[r]]
      e <- evidence[[m]][[r]] ## could use (height > 0) instead
      lapply(seq_along(unlist(o)), function(a)set.finding(domain, o[a], e[,a]))
    }

    out <- list()
    for (r in observed[[m]]){
      e <- evidence[[m]][[r]]
      out[[r]] <- do.call("rbind", lapply(seq_along(D[[r]]), function(a){
        retract(domain, O[[r]][a])            ## Remove conditioning on height[a]
        propagate(domain)                     ## Propagate evidence (on all, but O[a])

        ## MODIFIED FOR gRaven!
        ## Probabilities conditionally on heights[b], b!=a
        probs <- c(dat$allele[a], get.belief(domain, D[[r]][a]),
            get.belief(domain, Q[[r]][a]), r, dat[a, r+1])

        set.finding(domain, O[[r]][a], e[,a]) ## enter evidence for O[a], propagated under a+1
        probs

      }))
    }
    do.call("rbind", out)
  }

  f.prequential <- function(m){
    domain <- object$domains[[m]]
    D <- attr(domain, "D")
    O <- attr(domain, "O")
    Q <- attr(domain, "Q")
    dat <- object$data[[m]]

    ## No conditioning on peak heights (changes not yet propagated)
    ## retract(domain, unlist(attr(domain, "O")))

    f <- function(r, a){
      retract(domain, O[[r]][a])                             ## Remove evidence on O[a]
      propagate(domain)                                      ## propagate

      ## Modified for gRaven engine!
      ## Probabilities conditionally on heights[b], b!=a
      probs <- c(dat$allele[a], get.belief(domain, D[[r]][a]),
                 get.belief(domain, Q[[r]][a]), r, dat[a,r+1])

      set.finding(domain, O[[r]][a], evidence[[m]][[r]][,a]) ## re-enter evidence for O[a]

      probs
    }

    out <- matrix(NA, length(observed[[m]])*length(dat$allele), 7)
    i <- 0
    if (by.allele){ ## group by allele
      for (a in order(dat$allele)){
        for (r in observed[[m]]){
          i <- i+1
          out[i,] <- f(r,a)
        }
      }
    }
    else {
      for (r in observed[[m]]){ ## group by EPG
        for (a in order(dat$allele)){
          i <- i+1
          out[i,] <- f(r,a)
        }
      }
    }
    out
  }

  marker.probs <- switch(dist,
                         joint = f.joint,
                         conditional = f.conditional,
                         prequential = f.prequential,
                         known = f.known
                         )

  ## Return a list of data.frames
  out <- lapply(markers, function(m){
    d <- as.data.frame(marker.probs(m), stringsAsFactors = FALSE)
    names(d) <- c("allele", "seen", "unseen", "larger", "smaller", "trace", "height")
    d$marker <- m
    d
  })
  names(out) <- markers

  ## Note that evidence is left upon exit.
  out
}

