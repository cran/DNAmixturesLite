##' Loglikelihood function for DNA mixture analysis.
##'
##' @description
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which
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
##' The function \code{logL} is used to produce a log likelihood
##' function for one or more traces of DNA. \code{logL} calls one of
##' four other likelihood functions according to whether peak heights
##' or only peak presence/absence is used as observations, and also
##' whether there are any unknown contributors in the model.
##'
##' @note In the presence of unknown contributors, the returned
##' likelihood function relies on the Bayesian networks in the
##' \code{DNAmixture} object. If any evidence is entered or propagated
##' in these networks after creating the likelihood function, the
##' function will no longer compute the correct likelihood. In
##' particular, be ware of calling other functions in between calls to
##' the likelihood function, as they may change the state of the
##' networks.
##'
##'
##' @param mixture A \code{\link{DNAmixture}} model.
##' @param presence.only Set to TRUE to ignore peak heights and
##' evaluate the likelihood function considering peak presence and
##' absence (heights above and below threshold) only. Defaults to
##' FALSE.
##' @param initialize By default all entered
##' evidence is removed from the networks in \code{object}. Setting
##' \code{initialize = FALSE} should be done with care, and it is up
##' to the user to ensure that the likelihood being computed is meaningful.
##' @return A function, which takes a \code{\link{mixpar}} model parameter as argument.
##' @export
##' @author Therese Graversen
##' @examples
##' data(MC18, USCaucasian)
##' mixHp <- DNAmixture(list(MC18), k = 3, K = c("K1", "K2", "K3"), C = list(50),
##'                     database = USCaucasian)
##' p <- mixpar(rho = list(30), eta = list(34), xi = list(0.08),
##'             phi = list(c(K1 = 0.71, K3 = 0.1, K2 = 0.19)))
##' l <- logL(mixHp)
##' l(p)
logL <- function(mixture, presence.only = FALSE, initialize = TRUE){
  if (mixture$n.unknown > 0){
      if (presence.only) logLpres.UK(mixture, initialize = initialize)
      else logL.UK(mixture, initialize = initialize)
  }
  else {
      if (presence.only) logLpres.K(mixture)
      else logL.K(mixture)
  }
}

##' @export
##' @rdname logL
logL.UK <- function(mixture, initialize = TRUE){

  if (initialize) lapply(mixture$domains, initialize.domain)
  C <- mixture$C
  n.unknown <- mixture$n.unknown
  U <- mixture$U
  K <- mixture$K
  n_K <- lapply(mixture$data, function(d)subset(d, select = K))

  function(pararray){

    logL.m <- function(m){

      ## marker specific values
      domain <- mixture$domains[[m]]
      d <- mixture$data[[m]]
      rho <- pararray[,"rho"]
      xi <- pararray[,"xi"]
      eta <- pararray[,"eta"]
      phi <- pararray[,"phi"]

      for (r in seq_len(mixture$ntraces)){
        if (!all(is.na(d[,r+1]))){ ## if some heights are observed at this marker
          ## Set CPTs using parameter values
          evidence <- setCPT.O(domain, rho[[r]], xi[[r]], eta[[r]], phi[[r]][c(U, K)],
                               d[,r+1], C[[r]], n.unknown, n_K[[m]], attr(domain, "O")[[r]],
                               d$gets_stutter, d$can_stutter, d$stutter.from)
          ## Set evidence on O nodes
          lapply(seq_along(attr(domain, "O")[[r]]),
                 function(i)set.finding(domain, attr(domain, "O")[[r]][i], evidence[,i]))
        }
      }
      ## Propagate and get new log-normalizing constant
      if(is(domain, "RHuginDomain")) {propagate(domain)}
      get.normalization.constant(domain, log = TRUE)
    }

    ## Add up contributions from the M markers
    sum(sapply(seq_along(mixture$domains), logL.m))
  }
}


##' @export
##' @rdname logL
logL.K <- function(mixture){

  C <- mixture$C
  k <- mixture$k
  data <- mixture$data
  n_K <- lapply(data, function(d)subset(d, select = mixture$K))

  function(pars){
    rhos <- pars[,"rho"]
    etas <- pars[,"eta"]
    xis <- pars[,"xi"]
    phis <- pars[,"phi"]

    logL.m <- function(m){
      ## marker specific values
      d <- data[[m]]
      heights <- d[,seq_len(mixture$ntraces) + 1, drop = FALSE]
      can_stutter <- d$can_stutter
      gets_stutter <- d$gets_stutter
      stutter.from <- d$stutter.from
      n_K <- as.matrix(n_K[[m]])
      ##n_K_phi <- as.matrix(n_K[[m]]) %*% phi

      ## TODO: Swap the two loops
      one.allele <- function(a){

        one.trace <- function(rho, eta, xi, phi, height, threshold){
          ## Ensure that phi has the right ordering of contributors
          phi <- phi[mixture$K]

          ## As not all traces might have this allele:
          if (is.na(height)){return(0)}

          shape <- rho * (1 - xi*can_stutter[a]) * (n_K[a,] %*% phi)
          if (gets_stutter[a]){
            st <- stutter.from[a]
            ## can_stutter[st] = TRUE necessarily
            shape <- shape + rho * xi * (n_K[st,] %*% phi)
          }
          shape <- as.numeric(shape)

          ## likelihood for trace r, allele a
          if (height == 0)
            pgamma(threshold, shape = shape, scale = eta, log.p = TRUE)
          else
            dgamma(height, shape = shape, scale = eta, log = TRUE)
        }
        ## sum over traces, to get likelihood for allele a
        sum(mapply(FUN = one.trace, rhos, etas, xis, phis, heights[a,], C,
                   SIMPLIFY = TRUE))
      }
      ## sum over alleles, to get likelihood for marker m
      sum(sapply(seq_len(nrow(d)), one.allele))
    }
    ## sum over markers, to get overall likelihood
    sum(sapply(seq_along(mixture$markers), logL.m))
  }
}


##' @rdname logL
##' @export
logLpres.UK <- function(mixture, initialize = TRUE){

  if (initialize) lapply(mixture$domains, initialize.domain)

  C <- mixture$C
  n.unknown <- mixture$n.unknown
  U <- mixture$U
  K <- mixture$K
  n_K <- lapply(mixture$data, function(d)subset(d, select = K))

  function(pararray){

    logLpres.m <- function(m){

      ## marker specific values
      domain <- mixture$domains[[m]]
      d <- mixture$data[[m]]
      rho <- pararray[,"rho"]
      xi <- pararray[,"xi"]
      eta <- pararray[,"eta"]
      phi <- pararray[,"phi"]

      for (r in seq_len(mixture$ntraces)){
          if (!all(is.na(d[,r+1]))){ ## if some heights are observed at this marker
              ## Set CPTs using parameter values
              setCPT.D(domain, rho[[r]], xi[[r]], eta[[r]], phi[[r]][c(U, K)],
                        C[[r]], n.unknown, n_K[[m]], attr(domain, "D")[[r]],
                        d$gets_stutter, d$can_stutter, d$stutter.from)
              mapply(function(x, finding)set.finding(domain, x, finding), attr(domain, "D")[[r]], ifelse(d[,r+1] >= C[[r]], 0, 1))
          }
      }
      ## Propagate and get new log-normalizing constant
      if(is(domain, "RHuginDomain")) {propagate(domain)}
      get.normalization.constant(domain, log = TRUE)
    }

    ## Return the sum of contributions from the M markers
    sum(sapply(seq_along(mixture$domains), logLpres.m))
  }

}


##' @export
##' @rdname logL
logLpres.K <- function(mixture){

  domains <- mixture$domains
  C <- mixture$C
  k <- mixture$k
  data <- mixture$data
  n_K <- lapply(data, function(d)subset(d, select = mixture$K))

  function(pars){
    rhos <- pars[,"rho"]
    etas <- pars[,"eta"]
    xis <- pars[,"xi"]
    phis <- pars[,"phi"]

    logLpres.m <- function(m){
      ## marker specific values
      d <- data[[m]]
      heights <- d[,seq_len(mixture$ntraces) + 1, drop = FALSE]
      can_stutter <- d$can_stutter
      gets_stutter <- d$gets_stutter
      stutter.from <- d$stutter.from
      n_K <- as.matrix(n_K[[m]])

      one.allele <- function(a){

        one.trace <- function(rho, eta, xi, phi, height, threshold){
          ## Ensure that phi has the right ordering of contributors
          phi <- phi[mixture$K]

          ## As not all traces might have this allele:
          if (is.na(height)){return(0)}

          shape <- rho * (1 - xi*can_stutter[a]) * (n_K[a,] %*% phi)
          if (gets_stutter[a]){
            st <- stutter.from[a]
            ## can_stutter[st] = TRUE necessarily
            shape <- shape + rho * xi * (n_K[st,] %*% phi)
          }
          shape <- as.numeric(shape)

          ## likelihood for trace r, allele a
          if (height == 0)
            pgamma(threshold, shape = shape, scale = eta, log.p = TRUE, lower.tail = TRUE)
          else
            pgamma(threshold, shape = shape, scale = eta, log.p = TRUE, lower.tail = FALSE)
        }
        ## sum over traces, to get likelihood for allele a
        sum(mapply(FUN = one.trace, rhos, etas, xis, phis, heights[a,], C,
                   SIMPLIFY = TRUE))
      }
      ## sum over alleles, to get likelihood for marker m
      sum(sapply(seq_len(nrow(d)), one.allele))
    }
    ## sum over markers, to get overall likelihood
    sum(sapply(seq_along(mixture$markers), logLpres.m))
  }
}
