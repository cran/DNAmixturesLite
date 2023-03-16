##' Set tables for auxiliary variables \code{O}.
##'
##' @description Function for setting the tables on all nodes \code{O_r_a} for
##' mixture \code{r} in the network corresponding to one marker, using
##' peak heights and a set of parameters. The returned evidence can be
##' used for conditioning on observed peak heights.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##' \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @details The function sets probabilities in CPT according to
##' whether a peak is observed or not. If a peak is seen at allele
##' \eqn{a}, the conditional distribution is defined using the
##' p.d.f. \eqn{g} for the gamma distribution presented in
##' \code{\link{DNAmixtures}} as
##' \deqn{P(\texttt{O[a]} = \texttt{TRUE} | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = g(\texttt{heights[a]})/k_a
##' }{P(O[a] = TRUE | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = g(heights[a])/k_a}
##' and likelihood-evidence \eqn{(0, k_a)} for the states (\code{FALSE}, \code{TRUE}) is returned.
##' If the allele is not seen, the CPT is set using the c.d.f. \eqn{G} for the gamma distribution as
##' \deqn{P(\code{O[a]} = \code{FALSE} | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = G(\code{C}),
##' }{P(O[a] = FALSE | n_{ia}, n_{i,a+1}, i = 1,\ldots, k) = G(C),}
##' in which case likelihood-evidence (1,0) is returned for this allele.
##'
##' @section Technical comments: For the sake of speed and saving
##' memory, only probabilities for the CPT of \code{O[a]} are
##' calculated, and not the entire table.  The parent nodes are
##' \code{n_i_a} and possibly \code{n_i_(a+1)} for all unknown
##' contributors i; the \eqn{n_{ia}} for known contributors are
##' specified through the argument \code{n_K}.  The layout of the
##' table relies in particular on the fact that all nodes have a
##' statespace {0,1,2}, which this is *not* checked anywhere. Thus, if
##' other types of allele-count-nodes \code{n_i_a} are introduced
##' (such as could be relevant for Amelogenin), the function should be
##' changed accordingly.
##'
##' If invalid tables are set then any subsequent propagation will
##' fail. No roll-back functionality has so far been implemented to
##' fix this, and the easiest solution is to re-fit the mixture model.
##'
##' Note that the detection threshold and model parameters can be
##' specified as vectors, and so it is in theory possible to use allele-specific
##' values, if desired. This also holds for the functions
##' \code{\link{setCPT.D}} and \code{\link{setCPT.Q}}.

##' @export
##'
##' @param domain A \code{hugin.domain} modelling one marker
##' @param rho rho
##' @param xi xi
##' @param eta eta
##' @param phi phi
##' @param heights Peak heights
##' @param C Detection threshold
##' @param n.unknown Number of unknown contributors
##' @param n_K Allelecounts for known contributors
##' @param O Vector of names for binary nodes
##' @param gets_stutter Does the allele receive stutter? A boolean vector.
##' @param can_stutter Can the allele stutter?
##' @param stutter.from From which allele does the stutter come?
##' @return A matrix with each column containing the evidence to be set on \code{O[a]} when conditioning on the peak heights.
setCPT.O <- function(domain, rho, xi, eta, phi, heights, C, n.unknown,
                     n_K, O, gets_stutter, can_stutter, stutter.from){

  stopifnot(n.unknown > 0)
  ## Combinations of allele counts for all unknown contributors (at one allele)
  n_U_phi <- as.matrix(expand.grid(rep(list(0:2), n.unknown), KEEP.OUT.ATTRS = FALSE)) %*% head(phi, n.unknown)
  ## Allele counts for known contributors, all alleles
  n_K_phi <- if (!is.null(n_K)){as.matrix(n_K) %*% tail(phi, -n.unknown)} else rep(0, length(O))

  one.allele <- function(a){
    ## Shapes for allele a, contribution from allele a
    shapes <- rho * (1 - xi*can_stutter[a]) * (n_U_phi + n_K_phi[a])
    if (gets_stutter[a]){## Contribution from allele b
      b <- stutter.from[a]
      newshapes <- rho * xi * (n_U_phi + n_K_phi[b])
      shapes <- outer(shapes, newshapes, "+")
    }
    shapes <- as.vector(shapes)

    ## Compute the CPT frequencies.
    if (heights[a] > 0){
      ## If allele a is seen: P(O_a = TRUE) = f(height; allele counts)/Const
      ##      seen.densities <- dgamma(heights[a], shape = shapes, scale = eta)
      ##      Const <- max(seen.densities)
      ##      seen.probs <- seen.densities/Const

      ## Alternative
      log.seen.density <- dgamma(heights[a], shape = shapes, scale = eta, log = TRUE)
      logconst <- max(log.seen.density)
      seen.probs <- exp(log.seen.density - logconst)
      ## If Const is 0, corresponding to logconst = -Inf, then we should perhaps rather just stop.
      ## As it is now, the evidence will simply be constantly 0 and fail at propagation.
      Const <- exp(logconst)

      ## CPT alternates O_a = FALSE and O_a = TRUE
      cptfreqs <- as.numeric(rbind(1 - seen.probs, seen.probs))
      evidence <- c(0, Const) ## O_a is TRUE with weight Const
    }
    else {
      ## If height = 0: P(O_a = FALSE) = P(Z_a < C)
      unseen.probs <- ifelse(shapes > 0,
                             pgamma(C, shape = shapes, scale = eta),
                             1)
      cptfreqs <- as.numeric(rbind(unseen.probs, 1 - unseen.probs))
      evidence <- c(1, 0) ## O_a is FALSE
    }

    ## Set the CPT
    set.table(domain, O[a], cptfreqs, type = "cpt")

    ## Return evidence for future use in conditioning on peak heights
    evidence
  }

  ## Iterate over all alleles
  sapply(seq_along(O), one.allele)
}

