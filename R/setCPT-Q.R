##' Set tables for auxiliary variables \code{Q}
##'
##' @description Function for setting the tables on all nodes \code{Q_r_a} for
##' mixture \code{r} in the network corresponding to one marker, using
##' peak heights and a set of parameters.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##'
##' @details
##' The function sets conditional probabilities for auxiliary variables \code{Q_r_a},
##' which are used for finding probabilities of observing a smaller or
##' larger peak than the one observed for the DNA mixture.
##' The conditional probability tables are defined using the gamma c.d.f. as
##' \deqn{P(\texttt{Q[a]} = \texttt{TRUE} | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = G(\texttt{heights[a]}).
##' }{P(Q[a] = TRUE | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = G(heights[a]).}
##'
##' @export
##'
##' @param domain A \code{hugin.domain} modelling one marker
##' @param rho rho
##' @param xi xi
##' @param eta eta
##' @param phi phi ordered as \code{phi[c(U,K)]} where \code{U} and \code{K} are the names of unknown and known contributors.
##' @param heights Observed peak heights
##' @param n.unknown Number of unknown contributors
##' @param n_K Allelecounts for known contributors
##' @param Q Names for binary nodes
##' @param gets_stutter Does the allele receive stutter?
##' @param can_stutter Can the allele stutter?
##' @param stutter.from From which allele does the stutter come?
##' @return Scaling factors to be set as evidence when conditioning on the peak heights
##' @author Therese Graversen
##' @seealso For further details see \code{\link{setCPT.O}}.
setCPT.Q <- function(domain, rho, xi, eta, phi, heights, n.unknown,
                     n_K, Q, gets_stutter, can_stutter, stutter.from){

  stopifnot(n.unknown > 0)
  ## Combinations of allele counts for all unknown contributors (at one allele)
  n_U_phi <- as.matrix(expand.grid(rep(list(0:2), n.unknown), KEEP.OUT.ATTRS = FALSE)) %*% head(phi, n.unknown)
  ## Allele counts for known contributors, all alleles
  n_K_phi <- if (!is.null(n_K)){as.matrix(n_K) %*% tail(phi, -n.unknown)} else rep(0, length(Q))

  one.allele <- function(a){

    ## Contribution from allele a
    shapes <- rho * (1 - xi*can_stutter[a]) * (n_U_phi + n_K_phi[a])
    if (gets_stutter[a]){## Contribution from allele b
      b <- stutter.from[a]
      newshapes <- rho * xi * (n_U_phi + n_K_phi[b])
      shapes <- outer(shapes, newshapes, "+")
    }
    shapes <- as.vector(shapes)

    p.shorter <- pgamma(heights[a], shape = shapes, scale = eta)
    p.higher <- pgamma(heights[a], shape = shapes, scale = eta, lower.tail = FALSE)

    ## CPT alternates Q_a = FALSE and Q_a = TRUE between each row
    ## starting from Q_a = FALSE
    cptfreqs <- as.numeric(rbind(p.higher, p.shorter))
    set.table(domain, Q[a], cptfreqs, type = "cpt")
  }
  sapply(seq_along(Q), one.allele)
  invisible(NULL)
}
