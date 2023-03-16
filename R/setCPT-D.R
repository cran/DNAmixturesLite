##' Set tables for auxiliary variables \code{D}
##'
##' @description Function for setting the tables on all nodes \code{D_r_a} for
##' mixture \code{r} in the network corresponding to one marker, using
##' peak heights and a set of parameters.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @details
##' The function sets conditional probabilities for auxiliary variables \code{D_r_a},
##' which are used for obtaining probabilities of a peak
##' falling above resp. below the detection threshold under the current
##' model.
##' The conditional probability tables are defined using the gamma c.d.f. as
##' \deqn{P(\texttt{D[a]} = \texttt{TRUE} | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = G(\texttt{C}).
##' }{P(D[a] = TRUE | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = G(C).}
##' @export
##'
##' @param domain A \code{hugin.domain} modelling one marker
##' @param rho rho
##' @param xi xi
##' @param eta eta
##' @param phi phi
##' @param C Detection threshold
##' @param n.unknown Number of unknown contributors
##' @param n_K Allelecounts for known contributors
##' @param D Names for binary nodes
##' @param gets_stutter Does the allele receive stutter?
##' @param can_stutter Can the allele stutter?
##' @param stutter.from From which allele does the stutter come?
##' @return invisibly NULL
##' @seealso For further details see \code{\link{setCPT.O}}.
setCPT.D <- function(domain, rho, xi, eta, phi, C, n.unknown,
                     n_K, D, gets_stutter, can_stutter, stutter.from){

  ## Combinations of allele counts for all unknown contributors (at one allele)
  n_U_phi <- as.matrix(expand.grid(rep(list(0:2), n.unknown), KEEP.OUT.ATTRS = FALSE)) %*% head(phi, n.unknown)
  ## Allele counts for known contributors, all alleles
  n_K_phi <- if (!is.null(n_K)){as.matrix(n_K) %*% tail(phi, -n.unknown)} else rep(0, length(D))

  one.allele <- function(a){
    shapes <- rho * (1 - xi*can_stutter[a]) * (n_U_phi + n_K_phi[a])
    if (gets_stutter[a]){## Contribution from allele b
      b <- stutter.from[a]
      newshapes <- rho * xi * (n_U_phi + n_K_phi[b])
      shapes <- outer(shapes, newshapes, "+")
    }
    shapes <- as.vector(shapes)

    ## Find shape parameters
    ## Contribution from allele a
#    n_a <- as.matrix(merge(n_U, n_K[a,])) ## = n_U if n_K is NULL
#    shapes <- rho * (1 - xi*can_stutter[a]) * (n_a %*% phi)#

#    if (gets_stutter[a]){
#      ## Adding contributions from stutter from allele b
#      b <- stutter.from[a]
#      n_b <- as.matrix(merge(n_U, n_K[b,])) ## using n_Ua = n_Ub
#      newshapes <- rho * xi*can_stutter[b] * (n_b %*% phi)
#      shapes <- outer(shapes, newshapes, "+")
#    }
#    shapes <- as.vector(shapes)

    ## This works even when shape = 0.
    unseen.probs <- pgamma(C, shape = shapes, scale = eta)
    seen.probs <- pgamma(C, shape = shapes, scale = eta, lower.tail = FALSE)

    ## CPT alternates D_a = FALSE and D_a = TRUE between each row
    ## starting from D_a = FALSE
    cptfreqs <- as.numeric(rbind(seen.probs, unseen.probs))
    set.table(domain, D[a], cptfreqs, type = "cpt")
  }
  sapply(seq_along(D), one.allele)
  invisible(NULL)
}


