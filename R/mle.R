##' Maximisation of the likelihood for one or more mixed traces of DNA
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
##' @details The pre-specified constraints for the model parameter \code{pars} are for each mixture \code{r}
##' \itemize{
##' \item \code{0 <= pars[[r, "rho"]] < Inf}
##' \item \code{0 <= pars[[r, "eta"]] < Inf}
##' \item \code{0 <= pars[[r, "xi"]] <= 1}
##' \item \code{0 <= pars[[r, "phi"]][i] <= 1}
##' \item \code{sum(pars[[r, "phi"]]) == 1}.
##' \item If there are 2 or more unknown contributors, then the mixture
##' proportions for the unknown contributors are ordered decreasingly
##' with the first DNA mixture determining the order.
##' }
##' Further constraints can be specified by the user; for this see examples below.
##' @param mixture A \code{\link{DNAmixture}} object.
##' @param pars A \code{\link{mixpar}} parameter used as a starting value for the optimisation.
##' @param constraints Equality constraint function as function of an array of parameters.
##' @param phi.eq Should the mixture proportions be the same for all
##' traces? Defaults to FALSE.
##' @param val Vector of values to be satisfied for the equality constraints.
##' @param trace Print the evaluations of the likelihood-function during optimisation?
##' @param order.unknowns Should unknown contributors be ordered
##' according to decreasing contributions? Defaults to TRUE.
##' @param ... Further arguments to be passed on to \code{\link[Rsolnp]{solnp}}.
##' @return A list containing
##'    \item{mle}{The (suggested) MLE.}
##'    \item{lik}{The log of the likelihood (log e).}
##' as well as the output from the optimisation.
##' @author Therese Graversen
##' @export
##' @examples
##'
##' data(MC15, USCaucasian)
##' mix <- DNAmixture(list(MC15), C = list(50), k = 3, K = c("K1", "K2", "K3"), database = USCaucasian)
##' startpar <- mixpar(rho = list(24), eta = list(37), xi = list(0.08),
##'                    phi = list(c(K3 = 0.15, K1 = 0.8, K2 = 0.05)))
##' ml <- mixML(mix, startpar)
##' ml$mle
##'
##' \donttest{
##' data(MC15, USCaucasian)
##' mix <- DNAmixture(list(MC15), C = list(50), k = 3, K = "K1", database = USCaucasian)
##' startpar <- mixpar(rho = list(24), eta = list(37), xi = list(0.08),
##' phi = list(c(U1 = 0.05, K1 = 0.8, U2 = 0.15)))
##' ml <- mixML(mix, startpar)
##' ml$mle
##' }
##'
##' \donttest{
##' ## The following advanced example has a model with two DNA samples
##' ## and various parameter restrictions.
##' ## Be aware that the computation is rather demanding and takes
##' ## a long time to run with this lite-version of DNAmixtures.
##' ## With DNAmixtures (based on HUGIN), it takes only about one minute.
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = "K1", database = USCaucasian)
##' startpar <- mixpar(rho = list(24, 25), eta = list(37, 40), xi = list(0.08, 0.1),
##'                    phi = list(c(U1 = 0.05, K1 = 0.7, U2 = 0.25),
##'                                 c(K1 = 0.7, U2 = 0.1, U1 = 0.2)))
##' eqxis <- function(x){ diff(unlist(x[,"xi"])) }
##' ## Note that for these two traces, we do not expect phi to be equal.
##' ## Here we set stutter equal for all traces
##' ml.diff <- mixML(mix, startpar, eqxis, val = 0, phi.eq = FALSE)
##' ## Equal mixture proportions across traces
##' ml.eqphi <- mixML(mix, startpar, eqxis, val = 0, phi.eq = TRUE)
##' ## Likelihood ratio test for equal mixture proportions
##' pchisq(-2*(ml.eqphi$lik - ml.diff$lik), df = 1, lower.tail = FALSE)
##' }
##' @seealso \code{\link{DNAmixture}}
mixML <- function(mixture, pars, constraints = NULL, phi.eq = FALSE, val = NULL, trace = FALSE,
                  order.unknowns = TRUE, ...){
  R <- mixture$ntraces
  k <- mixture$k
  U <- mixture$U
  K <- mixture$K
  contr <- c(U, K)
  n.unknown <- mixture$n.unknown

  ## Return the list of (named) phis
  x2phi <- function(x){
    ## list with one phi per trace
    phi <- if (phi.eq){
      rep(list(tail(x, -3*R)), R)
    } else {
      split(tail(x, -3*R), rep(1:R, each = k))
    }

    ## Adding names of contributors to each phi
    mapply(function(th, contributor){names(th) <- contributor; th},
           phi, rep(list(contr), times=R), SIMPLIFY = FALSE)
  }

  ## List with phi_U for each trace
  x2phiU <- function(x){
    lapply(x2phi(x), function(th){head(th, n.unknown)})
  }

  ## vector --> parameter; i.e. array of lists: rows = traces, columns = parameters
  x2arr <- function(x){
    arr <- array(list(NULL), dim = c(R, 4), dimnames = list(NULL, c("rho", "eta", "xi", "phi")))
    arr[1:(3*R)] <- as.list(head(x, 3*R))
    arr[,"phi"] <- x2phi(x)
    arr
  }

  ## parameter --> vector
  parlist2x <- function(parlist){
    rex <- unlist(parlist[,1:3], use.names = FALSE)
    if (phi.eq){
      ## use phi for trace 1
      phi <- parlist[[1,4]][contr]
    } else {
      phi <- unlist(lapply(parlist[,4], function(x)x[contr]), use.names = FALSE)
    }
    c(rex, phi)
  }


  ## loglikelihood function
  logl <- logL(mixture)

  ## combined likelihood
  funvals <- numeric(0) ## only used if tracing.
  minus.loglikelihood <- function(x){
    xs <- x2arr(x)
    if (trace) print.mixpar(xs)
    val <- -logl(xs)
    if (trace) print(-val)
    val
  }

  ## Parameter boundaries
  lb <- rep(0, times = 3*R + ifelse(phi.eq, k, k*R))
  ub <- rep(c(Inf, 1), times = c(2*R, R + ifelse(phi.eq, k, k*R)))

  ## equality constraint returning a vector of sums for each phi
  if (phi.eq){
    phi.sum.constraint <- function(x){sum(tail(x, -3*R))}
    eqB <- 1
  }
  else {
    phi.sum.constraint <- function(x){sapply(x2phi(x), sum)}
    eqB <- rep(1, R)
  }

  if (!missing(constraints)){
    eqfun <- function(x){c(phi.sum.constraint(x), do.call(constraints, list(x2arr(x))))}
    eqB <- c(eqB, val)
  } else {
    eqfun <- phi.sum.constraint
  }

  ## Inequality constraint for traces with 2 or more unknowns.
  ## First phi in decreasing order.
  ## -- note that all traces have the same number of unknowns here.
  if (n.unknown > 1 & order.unknowns){
    phi.symmetry <- function(x){
      phiU <- x2phiU(x)[[1]] ## Irelevant to sort as [U]?
      diff(phiU, lag = 1)
    }
    n.diffs <- (n.unknown - 1) ##*R
    ineqLB <- rep(-1, n.diffs)
    ineqUB <- rep(0, n.diffs)
  }
  else phi.symmetry <- ineqLB <- ineqUB <- NULL

  ## Startvalue for the minimisation
  x0 <- parlist2x(pars)

  soln <- solnp(x0, fun = minus.loglikelihood,
                LB = lb, UB = ub,
                eqfun = eqfun, eqB = eqB,
                ineqfun = phi.symmetry,
                ineqLB = ineqLB, ineqUB = ineqUB,
                control=list(trace=0),
                ...
                )

  est <- x2arr(soln$pars)
  class(est) <- "mixpar"
  val <- -tail(soln$value, 1)
  out <- list(mle = est, lik = val,
              funvals = funvals,
              starting.point = x0,
              minimization.output = soln)
}
