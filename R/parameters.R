##' Parameters for DNA mixture models
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
##' @details The mixture parameter is a two-way array of lists;
##' columns correspond to the four model parameters \code{rho},
##' \code{eta}, \code{xi}, and \code{phi}, and rows correspond to the
##' mixtures included in the model.
##'
##' The print method is currently somewhat specialised, in that it
##' assumes that \code{rho}, \code{eta}, and \code{xi} are merely real
##' numbers. \code{phi} is assumed to consist of a named vector per
##' mixture; the names, or order of names, can differ between mixtures.
##'
##' @author Therese Graversen
##' @export
##' @param rho Amplification factor.
##' @param eta Scale parameter in the gamma distribution.
##' @param xi Stutter parameter.
##' @param phi Named vector of the fraction of DNA contributed by each contributor.
##' @param parlist A list of parameters of class \code{mixpar}
##' @return An object of class "mixpar".
##' @examples
##' ## A parameter for two mixtures
##' q <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08, 0.08),
##'             phi = list(c(Anna = 0.5, Peter = 0.2, U1 = 0.3),
##'                        c(U1 = 0.5, Anna = 0.2, Peter = 0.3)))
##' ## Equivalent to specifying the parameter for each mixture and then combining.
##' p1 <- mixpar(rho = list(30), eta = list(30), xi = list(0.08),
##'              phi = list(c(Anna = 0.5, Peter = 0.2, U1 = 0.3)))
##' p2 <- mixpar(rho = list(30), eta = list(30), xi = list(0.08),
##'              phi = list(c(U1 = 0.5, Anna = 0.2, Peter = 0.3)))
##' p <- mixpar(parlist = list(p1, p2))
mixpar <- function(rho = NULL, eta = NULL, xi = NULL, phi = NULL, parlist = NULL){

  if (!missing(parlist)){
    stopifnot(all(sapply(parlist, function(x)class(x)=="mixpar")))
    rho <- lapply(parlist, function(x)x[, "rho"])
    eta <- lapply(parlist, function(x)x[, "eta"])
    xi <- lapply(parlist, function(x)x[, "xi"])
    phi <- lapply(parlist, function(x)x[, "phi"])
  }

  arr <- array(list(), dim = c(length(rho), 4),
               dimnames = list(NULL, c("rho", "eta", "xi", "phi")))

  arr[,1] <- rho
  arr[,2] <- eta
  arr[,3] <- xi
  arr[,4] <- phi
  structure(arr,
            class = "mixpar")
}

##' @rdname mixpar
##' @method print mixpar
##' @export
##' @param x An object of class \code{mixpar}.
##' @param scientific Should scientific notation be used?
##' @param digits The number of digits to print
##' @param ... arguments passed to print
print.mixpar <- function(x, digits = max(3L, getOption("digits") - 3L),
                         scientific = FALSE, ...){
  phi.names <- unique(do.call("c", lapply(x[,"phi"], names)))
  xx <- x
  xx[,"phi"] <- lapply(xx[,"phi"], function(v){out <- v[phi.names]; names(out) <- phi.names;out})
  print(format(as.data.frame(t(apply(xx, 1, unlist))), digits = digits, scientific = scientific),
        quote = FALSE, print.gap = 3L, ...)
  invisible(x)
}
