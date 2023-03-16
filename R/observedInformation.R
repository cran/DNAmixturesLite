##' Estimated asymptotic variance matrix for MLE
##'
##' @description
##' Provided that the user specifies the MLE as well as the
##' constraints used in the maximisation, this function computes an
##' estimate of the variance of the MLE based on the observed
##' information. The observed information is obtained by numerically
##' deriving the Hessian of the log-likelihood function.
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
##' @details As the user can apply highly customized constraints to the model
##' parameters when maximising with \code{\link{mixML}}, it is a
##' complicated matter to write a generic function for computing the
##' asymptotic variance matrix. We have thus restricted attention to
##' the case where each of the (multi-dimensional) parameters \code{rho},
##' \code{eta}, \code{xi} and \code{phi} can be either
##'
##' \itemize{
##' \item fixed at known values
##' \item unknown, but common across traces
##' \item unconstrained
##' }
##'
##' @param mixture A \code{\link{DNAmixture}}.
##' @param mle A \code{mixpar}, typically obtained by \code{mixML}.
##' @param npars A list of integers specifying the number of each of the four
##' parameters \code{rho}, \code{eta}, \code{xi} and \code{phi}. Allowed values are
##' \describe{
##' \item{\code{0}}{The parameter is fixed, but might differ across mixtures.}
##' \item{\code{1}}{The parameter is equal across mixtures.}
##' \item{\code{N}}{There is one parameter for each of the \code{N} mixtures in the model.}
##' }
##' @param method Method for numeric differentiation used in \code{\link[numDeriv]{hessian}}
##' @param ... Arguments to be passed on to other methods.
##' @return The mle and the estimated covariance matrix in different parametrisations
##' \item{\code{cov} and \code{mle}}{\eqn{\rho, \eta, \xi, \phi}}
##' \item{\code{cov.trans} and \code{mle.trans}}{\eqn{\mu, \sigma, \xi, \phi}}
##' \item{\code{cov.res}}{A non-singular covariance matrix for a
##' reparametrisation of \eqn{\rho, \eta, \xi, \phi}, collapsing the
##' parameters according to the specified constrains and removing one
##' contributor from \eqn{\phi}.}
##'
##' An integer suffix is used to indicate which mixture the parameter
##' is associated with. In the restricted covariance matrix, all fixed
##' parameters are left out. If parameters are equal accross mixtures,
##' the suffix for this parameter will be \code{.1}. If parameters are
##' unconstrained and there are \code{N} mixtures in the model,
##' suffixes are \code{.1, \ldots, .N}
##' @export
##' @examples
##'
##' data(MC18, USCaucasian)
##' mixHp <- DNAmixture(list(MC18), k = 3, K = c("K1", "K2", "K3"), C = list(50),
##'                     database = USCaucasian)
##' p <- mixpar(rho = list(30), eta = list(34), xi = list(0.08),
##'             phi = list(c(K1 = 0.71, K3 = 0.1, K2 = 0.19)))
##' mlHp <- mixML(mixHp, pars = p)
##' ## Find the estimated covariance matrix of the MLE
##' V.Hp <- varEst(mixHp, mlHp$mle, npars = list(rho=1,eta=1,xi=1,phi=1))
##' V.Hp$cov ## using (rho, eta)
##' V.Hp$cov.trans ## using (mu, sigma)
##' ## The summary is a table containing the MLE and their standard errors
##' summary(V.Hp)
##'
##' \donttest{
##' data(MC18, USCaucasian)
##' mixmult <- DNAmixture(list(MC18), C = list(50), k = 3, K = c("K1", "K2"), database = USCaucasian)
##' startpar <- mixpar(rho = list(30), eta = list(28), xi = list(0.08),
##'                    phi = list(c(U1 = 0.2, K1 = 0.7, K2 = 0.1)))
##' ml.mult <- mixML(mixmult, startpar)
##' Vmult <- varEst(mixmult, ml.mult$mle, list(rho=1,eta=1,xi=1,phi=1))
##' summary(Vmult)
##' }
##'
##' \donttest{
##' ## Be aware that the following two advanced examples are computationally demanding and
##' ## typically have a runtime of several minutes with the lite-version of DNAmixtures.
##'
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = c("K1", "K2"),
##' database = USCaucasian)
##' startpar <- mixpar(rho = list(30, 30), eta = list(28, 35), xi = list(0.08, 0.1),
##'                    phi = list(c(U1 = 0.05, K1 = 0.7, K2 = 0.25),
##'                                 c(K1 = 0.7, K2 = 0.1, U1 = 0.2)))
##' eqxis <- function(x){ diff(unlist(x[,"xi"])) }
##' ## Here we set stutter equal for all traces
##' ml.diff <- mixML(mix, startpar, eqxis, val = 0, phi.eq = FALSE)
##' V.diff <- varEst(mix, ml.diff$mle, list(rho=2,eta=2,xi=1,phi=2))
##' summary(V.diff)
##'
##' ## Fixing stutter to 0.07
##' xival <- function(x){unlist(x[,"xi"])}
##' ml.eq <- mixML(mix, startpar, xival, val = c(0.07, 0.07), phi.eq = FALSE)
##' V.eq <- varEst(mix, ml.eq$mle, list(rho=2,eta=2,xi=0,phi=2))
##' summary(V.eq)
##' }
varEst <- function(mixture, mle, npars, method = "Richardson", ...){
  R <- mixture$ntraces
  k <- mixture$k
  mle[,4] <- lapply(mle[,4], function(x)x[c(mixture$U, mixture$K)])
  npars <- unlist(npars[dimnames(mle)[[2]]])
  if (!all(npars %in% c(0, 1, R))){stop(paste("The number of each parameter has to be", 0, ",", "1", "or", R))}

  fixed <- which(npars == 0)
  equal <- which(npars == 1)
  npars[[4]] <- npars[[4]]*(mixture$k-1)
  xlen <- sum(npars)

  extend.phi <- function(y){
    y <- c(y, 1-sum(y))
    names(y) <- c(mixture$U, mixture$K)
    y
  }

  logit <- function(v)log(v)-log(1-v)
  expit <- function(v)exp(v)/(1+exp(v))
  f <- function(i){
    switch(i,
           log,
           log,
           logit,
           function(th)log(th)-log(1-sum(th))
           )
  }

  f.inv <- function(i){
    switch(i,
           exp,
           exp,
           expit,
           function(y) exp(y)/(1 + sum(exp(y)))
           )
  }

  x2pars <- function(x){
    xs <- x
    p <- array(list(NA), dim = c(R, 4), dimnames = dimnames(mle))
    for (i in 1:4){

      if (npars[i] > 0){ ## Parameter i is not fixed
        ## i'th parameters
        v <- head(xs, npars[i])
        ## Rest of untransformed parameters
        xs <- tail(xs, -npars[i])

        ## Equal parameters, so all ready, except phi that needs extension
        if (i %in% equal){
          ## Transform to original parametrisation
          v <- f.inv(i)(v)
          if (i == 4){
            v <- extend.phi(v)
          }
          ## duplicate the parameter to all traces
          p[,i] <- list(v)
        }
        else {
          ## Different parameters
          if (i == 4){
            ## R chunks of size k-1
            p[,i] <- split(v, rep(1:(npars[i]/(k-1)), each = k-1))
            ## Extend to size k
            p[,i] <- lapply(p[,i], function(y)extend.phi(f.inv(4)(y)))
          }
          else{
            ## R chunks, each of size 1
            p[,i] <- split(v, 1:npars[i])
            p[,i] <- lapply(p[,i], f.inv(i))
          }
        }
      }
      else { ## fixed, so use mle directly
        p[,i] <- mle[,i]
      }
    }
    p
  }

  pars2x <- function(p){
    x <- list(4)
    for (i in 1:4){
      if (!(i %in% fixed)){
        if (i %in% equal){
          x[[i]] <- p[1,i] ## one parameter
        }
        else {
          x[[i]] <- p[,i] ## R parameters
        }
        if (i == 4) x[[4]] <- lapply(x[[4]], head, -1)
        x[[i]] <- lapply(x[[i]], f(i))
      }
    }
    unlist(x)
  }

  ## Gradient of the transformation in the MLE
  ## The phi matrix becomes
  ## phi*I - {phi_i * phi_j}_ij
  ## And (rho, eta, xi) are coordinatewise transformations

  mle.trans <- pars2x(mle)
  xs <- mle.trans
  g <- rep(list(numeric(0)), 4)

  for (i in 1:3){
    if (npars[i] > 0){
      g[[i]] <- head(xs, npars[i])
      xs <- tail(xs, -npars[i])
    }
  }
  if (npars[4] > 0){
    g[[4]] <- split(xs, rep(1:(npars[4]/(k-1)), each = (k-1)))
    ths <- lapply(g[[4]], function(y) exp(y)/(1 + sum(exp(y))))
  }

  ## This should work even if some matrices are 0x0
  elts <- c(exp(g[[1]]), exp(g[[2]]), expit(g[[3]])-expit(g[[3]])^2)
  D1 <- diag(elts, nrow = length(elts))
  Ds <- lapply(ths, function(y)diag(y, nrow = length(y))-matrix(y, ncol=1)%*%matrix(y, nrow=1))
  G <- bdiag(c(list(D1), Ds))

  logL <- logL(mixture)
  l <- function(x){
    -logL(x2pars(x))
  }

  ## Dimension names
  dns.res <- mapply(paste,
                    c(rep(c("rho", "eta", "xi"), times = npars[1:3]), rep(head(names(mle[[1,4]]),-1), times = npars[4]/(k-1))),
                    c(seq_len(npars[1]), seq_len(npars[2]), seq_len(npars[3]), rep(seq_len(npars[4]/(k-1)), each = k-1)),
                    MoreArgs = list(sep = "."))

  last <- tail(names(mle[[1,4]]),1) ## name of last contributor
  dns.trans <- mapply(paste,
                      c(rep(c("log(rho)", "log(eta)", "logit(xi)"), times = npars[1:3]),
                        paste0("log(", rep(head(names(mle[[1,4]]),-1), times = npars[4]/(k-1)),"/",last,")")),
                      c(seq_len(npars[1]), seq_len(npars[2]), seq_len(npars[3]), rep(seq_len(npars[4]/(k-1)), each = k-1)),
                      MoreArgs = list(sep = "."))

  ## The hessian of the transformed problem
  H <- hessian(l, pars2x(mle), method = method)
  dimnames(H) <- dimnames(G) <- list(dns.trans, dns.trans)

  ## Covariance matrix in normal parametrisation, with dimension names
  cov.res <- as.matrix(G %*% solve(H) %*% G)
  dimnames(cov.res) <- list(dns.res, dns.res)

  ## Extend to full parametrisation with a set of parameters per trace
  ## Also, here phi is extended to include all contributors
  A <- matrix(0, nrow = nrow(cov.res), ncol = (3+k)*R)
  i <- 1
  for (j in 1:3L){
    if (npars[j] == 1){
      A[i, ((j-1)*R+1):(j*R)] <- matrix(1, nrow = 1, ncol = R)
    }
    else if (npars[j] == R){
      A[i:(i+R-1), ((j-1)*R + 1):(j*R)] <- diag(1, nrow = R, ncol = R)
    }
    i <- i+npars[j]
   }
  if (npars[4] == 1*(k-1)){
    A[i:nrow(A), (3*R+1):ncol(A)] <- do.call(cbind, rep(list(cbind(diag(k-1),matrix(-1,nrow=k-1,ncol=1))), R))
  }
  else if (npars[4] == R*(mixture$k-1)){
    A[i:nrow(A), (3*R+1):ncol(A)] <- as.matrix(bdiag(rep(list(cbind(diag(k-1),matrix(-1,nrow=k-1,ncol=1))), R)))
  }

  ## A <- NULL
  cov <- t(A) %*% cov.res %*% A
  dns <- paste(unlist(rep(list("rho", "eta", "xi", paste0("phi.", names(mle[[1,4]]))), each = R)),
               c(rep(seq_len(R), times=3), rep(seq_len(R), each = k)),
               sep = ".")
  dimnames(cov) <- list(dns, dns)

  ## Change to (mu, sigma)-parametrisation
  B <- diag(1, nrow = nrow(cov), ncol = ncol(cov))
  for (i in 1:R) {## Rho row
    B[i,i] <- mle[[i,"eta"]] ## Mu column
    B[i,i+R] <- -0.5*mle[[i,"rho"]]^(-3/2) ## Sigma column
  }
  for (i in (R+1):(2*R)){ ## Eta row
      B[i,i-R] <- mle[[i-R,"rho"]] ## Mu column
      B[i,i] <- 0 ## Sigma column
  }
  cov.trans <- t(B) %*% cov %*% B
  dns.trans <- paste(unlist(rep(list("mu", "sigma", "xi", paste0("phi.", names(mle[[1,4]]))), each = R)),
                     c(rep(seq_len(R), times=3), rep(seq_len(R), each = k)),
                     sep = ".")
  dimnames(cov.trans) <- list(dns.trans, dns.trans)
  dimnames(A) <- list(dns.res, dns)
  dimnames(B) <- list(dns, dns.trans)

  mle.trans <- mle
  dimnames(mle.trans)[[2]][1:2] <- c("mu", "sigma")
  mle.trans[,1] <- apply(mle, 1, function(x)x[[1]]*x[[2]])
  mle.trans[,2] <- lapply(mle[,1], function(x)1/sqrt(x))

  structure(
    list(cov = cov, mle = mle, cov.res = cov.res, mle.trans = mle.trans, cov.trans = cov.trans, hessian.raw = H),
    class = "mixVarEst", A = A, B = B, f = f, f.inv = f.inv, G = G
    )
}

##' @export
##' @method summary mixVarEst
##' @param object An object of class \code{mixVarEst}, typically obtained by a call to \code{\link{varEst}}.
##' @param transform Should the parameterisation \eqn{(\mu, \sigma)} be used? Defaults to \code{FALSE}.
##' @rdname varEst
summary.mixVarEst <- function(object, transform = FALSE, ...){
  if (transform){
    est <- object$mle.trans
    SE <- sqrt(diag(object$cov.trans))
  }
  else {
    est <- object$mle
    SE <- sqrt(diag(object$cov))
  }
  structure(data.frame(Estimate = unlist(est),
                       StdErr = SE),
            class = c("summary.mixVarEst", "data.frame"))
}

##' @rdname varEst
##' @export
##' @method print summary.mixVarEst
##' @param x An object of class \code{"summary.mixVarEst"}.
##' @param digits Number of significant digits to print
##' @param print.gap Distance between columns in the printing of the summary.
##' @param scientific Should scientific notation be used?
print.summary.mixVarEst <- function(x, digits = max(3, getOption("digits") - 3), scientific = FALSE, print.gap = 3L, ...){
  print(format(x, digits = digits, scientific = FALSE, ...), quote = FALSE, print.gap = print.gap)
  invisible(x)
}
