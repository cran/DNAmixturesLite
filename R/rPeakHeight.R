##' Simulate peak heights from a DNA mixture model.
##'
##' @description A set of peak heights is sampled for each mixture in the
##' model. Note that if the mixtures cover different sets of markers,
##' so will the simulated peak heights. If there are unknown
##' contributors, their genotypes are sampled from the networks in
##' \command{mixture$domains}, but these genotypes are not stored.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @param mixture A \code{DNAmixture}.
##' @param pars A \code{mixpar} model parameter.
##' @param nsim Number of simulations for each allele.
##' @param markers Subset of markers to simulate
##' @param dist How should observed peak heights be taken into account?
##' Possible choices are
##' \describe{
##' \item{\code{"joint"}}{ which does not take peak heights into account.}
##' \item{\code{"conditional"}}{ which takes all peak heights into account, except that for the peak under consideration.}
##' }
##' @param use.threshold Should peak heights under the detection
##' threshold C be set to 0? Defaults to \code{TRUE}, corresponding to simulating under the model fitted.
##' @param initialize Should all propagated evidence be removed from
##' the network? Defaults to \code{TRUE}. If set to \code{FALSE}, the user should
##' make sure to check whether the networks represent the distribution
##' of interest.
##' @param ... Not used.
##' @return A list of one three-way array per marker; the three margins correspond to traces, alleles, and simulations.
##' @author Therese Graversen
##' @examples
##' data(MC15, MC18, USCaucasian)
##'
##' mixP <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = c("K1", "K3", "K2"),
##'                    database = USCaucasian)
##' pP <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08, 0.08),
##'              phi = list(c(K2 = 0.1, K3 = 0.2, K1 = 0.7), c(K2 = 0.1, K3 = 0.2, K1 = 0.7)))
##' rpk <- rPeakHeight(mixP, pP, nsim = 5)
##'
##' \donttest{
##' mixD <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = c("K1", "K3"),
##'                    database = USCaucasian)
##' pD <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08, 0.08),
##'              phi = list(c(U1 = 0.1, K3 = 0.2, K1 = 0.7), c(U1 = 0.1, K3 = 0.2, K1 = 0.7)))
##' rpj <- rPeakHeight(mixD, pD, nsim = 5, dist = "joint")
##' rpc <- rPeakHeight(mixD, pD, nsim = 5, dist = "conditional")
##' }
##' @export
rPeakHeight <- function(mixture, pars, nsim = 1, markers = mixture$markers,
                        dist = c("joint", "conditional"),
                        use.threshold = TRUE, initialize = TRUE, ...){
  rhos <- pars[,"rho"]
  etas <- pars[,"eta"]
  xis <- pars[,"xi"]
  phis <- pars[,"phi"]
  phi_U <- lapply(phis, function(t)t[mixture$U])
  phi_K <- lapply(phis, function(t)t[mixture$K])

  n.unknown <- mixture$n.unknown
  if (n.unknown == 0) dist <- "known"
  else {
    dist <- match.arg(dist)
    ## Remove any propagated evidence
    if (initialize) lapply(mixture$domains, initialize.domain)
  }

  ## simulate from network, exctract allele counts, convert to numeric
  ## values (from factors) Note that in many cases we only need
  ## to simulate for two alleles, not the entire network.
  rAlleleCount <- function(d, nodes){
    ns <- subset(simulate(d, nsim = nsim), select = nodes)
    as.data.frame(lapply(ns, function(f)as.numeric(levels(f))[f]))
  }

  ## Shape contributions from known contributors
  if (length(mixture$K) > 0){
    shapes_K <- getShapes(mixture, pars)
  }

  f.known <- function(m) {
    dim(shapes_K[[m]]) <- c(dim(shapes_K[[m]]), 1)
    shapes_K[[m]]
  }

  f.joint <- function(m){
    domain <- mixture$domains[[m]]
    d <- mixture$data[[m]]

    ## Simulate allele counts for unknown contributors
    ac <- as.matrix(rAlleleCount(domain, nodes = attr(domain, "n")))

    one.trace <- function(rho, xi, tU){
      ## |alleles| x |alleles| matrix
      S <- diag(1 - xi*d$can_stutter, nrow = length(d$can_stutter))
      S <- S + xi*(col(S) == d$stutter.from[row(S)])
      do.call("cbind", ## one column per simulation
              lapply(1:nsim,
                     function(i)rho*(S %*% matrix(ac[i,], ncol = n.unknown) %*% tU)))
    }
    shapes <- array(dim = c(mixture$ntraces, nrow(d), nsim))
    for (r in mixture$observed[[m]]){
      shapes[r,,] <- one.trace(rhos[[r]], xis[[r]], phi_U[[r]])
    }
    if (length(mixture$K) > 0){
      ## Add shape_K () to each simulated set of shapes for each value
      ## of margin 3.
      shapes <- sweep(shapes, MARGIN = c(1,2), shapes_K[[m]], '+')
    }
    shapes
  }

  f.conditional <- function(m){
    domain <- mixture$domains[[m]]
    d <- mixture$data[[m]]
    O <- attr(domain, "O")

    ## Set cpts for the mixture
    evidence <- setCPT(mixture, pars, markers = m)[[1]]
    for (r in mixture$observed[[m]]){
      for (a in seq_along(O)){
        set.finding(mixture$domain[[m]], O[[r]][a], evidence[[r]][,a])
      }
    }

    ## Shapes array initialising with NA's
    shapes <- array(dim = c(mixture$ntraces, nrow(d), nsim))

    ## For each trace that includes marker m fill in shapes.
    for (r in mixture$observed[[m]]){

      ## Run along alleles of marker m one by one
      shapes[r,,] <- do.call(rbind, lapply(seq_along(O[[r]]), function(a){
        retract(domain, O[[r]][a]) ## Remove conditioning on height[a]
        propagate(domain)

        ## Get the shapes for unknown contributors (same U for all traces)
        na <- attr(domain, "n")[a,]
        if (d$gets_stutter[a]){
          nb <- attr(domain, "n")[d$stutter.from[a],]
          ac <- rAlleleCount(domain, nodes = c(na, nb))
          n_U_a <- as.matrix(subset(ac, select = na))
          n_U_b <- as.matrix(subset(ac, select = nb))
          s <- rhos[[r]] * ( (1 - xis[[r]]*d$can_stutter[a]) * (n_U_a %*% phi_U[[r]]) + xis[[r]] * (n_U_b %*% phi_U[[r]]))
        }
        else {
          n_U_a <- as.matrix(rAlleleCount(domain, nodes = na))
          s <- rhos[[r]] * (1 - xis[[r]]*d$can_stutter[a]) * (n_U_a %*% phi_U[[r]])
        }
        ## re-enter observation a (not propagated, in effect next iteration)
        set.finding(domain, O[[r]][a], evidence[[r]][,a])

        dim(s) <- NULL
        s
      }))
    }
    if (length(mixture$K) > 0){
      ## Add shape_K () to each simulated set of shapes for each value
      ## of margin 3.
      shapes <- sweep(shapes, MARGIN = c(1,2), shapes_K[[m]], '+')
    }

    shapes
  }

  ## Choose appropriate function for calculating simulated shapes
  shapeFun <- switch(dist,
                     joint = f.joint,
                     conditional = f.conditional,
                     known = f.known
                     )

  ## Return a list of arrays with dimension |traces| x |alleles| x nsim
  out <- lapply(markers, function(m){
    A <- nrow(mixture$data[[m]])
    arr <- array(dim = c(mixture$ntraces, A, nsim))
    sh <- shapeFun(m)

    for (r in mixture$observed[[m]]){
      for (a in seq_len(A)){
        arr[r, a, ] <- rgamma(nsim, shape = sh[r,a,], scale = etas[[r]])
        arr[r, a, ] <- ifelse(arr[r, a, ] < mixture$C[[r]], 0, arr[r, a, ])
      }
    }
    arr
  })
  names(out) <- markers

  ## Note that propagated evidence is left on exit.
  out
}


