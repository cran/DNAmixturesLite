##' Create a DNA mixture model
##'
##' @description
##' A model object for analysis of one or more DNA mixtures. For a
##' brief overview of the package functionality, see
##' in particular \code{\link{DNAmixtures}}.
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
##' @details
##'
##' The names for known contributors can be chosen freely, whereas
##' unknown contributors are always termed \code{U1,U2, ...}.
##'
##' We allow for stutter to an allele one repeat number shorter. The
##' range of alleles at a marker is defined by the union of alleles
##' specified though the peak heights, the reference profiles, and the
##' allele frequencies. Any alleles that are included, but not found
##' in the database, will be assigned frequency \code{NA}, and it is then up
##' to the user to decide on further actions. If a particular mixture
##' has no observations at a marker, the heights are set to \code{NA}, but if
##' the mixture has some peaks at that marker, then missing heights
##' are all set to 0. Note that we hereby cover the possibility that
##' mixtures are analysed with different kits, and so are observed at
##' different markers. We do not (readily) allow kits to have
##' different ranges of possible alleles at one marker.
##'
##' If amelogenin is included in the analysis, the marker should be
##' named \code{"AMEL"} and an integer coding such as X=0, Y=1, where X is
##' assigned a lower number than Y, should be used. Note that in terms
##' of amelogenin, the allele frequencies have a slighly different
##' interpretation to that for other markers, in that they denote the
##' probability of having an \emph{additional} X or Y to the X that
##' all people have. Thus, a natural choice will be \eqn{p(X)=p(Y)=0.5},
##' denoting equal probability of a male or female contributor.
##'
##' @param data A list containing one \code{data.frame} for each DNA
##' mixture. Note, that in the special case of analysing just one
##' mixture, this still has to be specified as list(data). Each
##' dataset should contain variables \code{marker}, \code{allele}, and
##' \code{frequency}. Optionally, also a column for each reference
##' profile specified in \code{K}.

##' @param k Number of contributors.
##' @param C A list of thresholds, one for each mixture.
##' @param K Names of reference profiles; these can be chosen freely,
##' but should match (possibly only a subset of) the names specified
##' by the reference profiles.
##' @param database A data.frame containing at least variables \code{marker}, \code{allele}, \code{frequency}.
##' @param reference.profiles A data.frame containing allele counts for each reference profile, if not specified in \code{data}.
##' @param dir Location of network files if loading or saving the networks.
##' @param domainnamelist Names of marker-wise network files (without hkb-extension). Default is the set of markers.
##' @param load Read networks from disk?
##' @param write Save networks as hkb files?
##' @param dyes A list containing a list of dyes indexed by markers
##' @param triangulate Triangulate the networks? Default is to
##' triangulate the network using a good elimination order.
##' @param compile Compile the networks?
##' @param compress Compress the network? Defaults to \code{TRUE} and is
##' strongly recommended for models with a large number of
##' contributors.
##' @param use.order Should the default elimination order be used for triangulation?
##' Otherwise the "total.weight" heuristic for triangulation in Hugin is used.
##' @return An object of class \code{DNAmixture}. This contains amongst other things
##' \item{markers}{The joint set of markers used for the mixtures specified.}
##' \item{domains}{For models involving unknown contributors,
##' a list containing one Bayesian network (\code{hugin.domain}) per marker;
##' see \code{\link{buildMixtureDomains}} for details on the networks}
##' \item{data}{A list containing for each marker the combined allele frequencies,
##' peak heights, and reference profiles as produced by \code{\link{DNAmixtureData}}.}
##' @export
##' @examples
##' data(MC15, MC18, USCaucasian)
##' DNAmixture(list(MC15, MC18), C = list(50,50), k = 3, K = c("K1", "K2"),
##'            database = USCaucasian)
##' DNAmixture(list(MC15, MC18), C = list(50,50), k = 3, K = c("K3", "K1", "K2"),
##'            database = USCaucasian)
DNAmixture <- function(data, k, C, database, K = character(0),
                       reference.profiles = NULL,
                       dir = character(0), domainnamelist = NULL,
                       load = FALSE, write = FALSE,
                       dyes = NULL, triangulate = TRUE, compile = TRUE, compress = TRUE, use.order = TRUE){

  ## Names of the datasets used, to be stores as the "call"
  call <- deparse(substitute(data))

  n.unknown <- k - length(K)
  if (n.unknown < 0) stop("Too many known contributors specified")
  if (length(C) != length(data)) stop("There needs to be one threshold per trace")
  if (!is.null(dyes) & length(dyes) != length(data)) stop("Specify a dye per trace")
  stopifnot(class(data)=="list")

  ntraces <- length(data)

  ## Any peaks being censored by the chosen threshold?
  for (r in 1:ntraces){
    if (any((data[[r]]$height < C[[r]]) & (data[[r]]$height > 0)))
      warning("Some observed peaks are below the given detection threshold and are set to 0",
              call. = FALSE, immediate. = TRUE)
    data[[r]]$height <- ifelse(data[[r]]$height < C[[r]], 0, data[[r]]$height)
  }

  data <- DNAmixtureData(data, database = database, reference.profiles = reference.profiles, K = K)
  observed <- lapply(data, function(d)as.vector(which(apply(d[,1:ntraces + 1, drop=F], 2, function(x)any(!is.na(x))))))
#  lapply(data, function(d)as.vector(which(sapply(d[,1:ntraces + 1], function(x)any(!is.na(x))))))

  out <- list(
    ## Markers
    markers = names(data),

    ## Data
    data = data,

    ## Number of traces
    ntraces = ntraces,

    ## Number of contributors
    k = k,

    ## Number of unknown contributors
    n.unknown = n.unknown,

    ## Names of known contributors
    K = K,

    ## Detection threshold.
    C = C,

    ## Which traces are observed for given markers
    observed = observed,

    ## Dyes used in the kit
    dyes = dyes,

    ## names of the mixtures included
    call = call
    )

  ## If there are unknown contributors we need networks
  if (n.unknown > 0){

    ## Names for any unknown contributors. Used in labeling parameters
    out$U <- paste("U", seq_len(n.unknown), sep = "")

    ## Default directory is of the form xU/ for x unknown contributors
    if (length(dir) == 0){
      dir <- paste(n.unknown, "U/", sep = "")
    }

    ## Default domain names are the marker names
    if (is.null(domainnamelist)){
      domainnamelist <- as.list(names(data))
      names(domainnamelist) <- names(data)
    }

    if (!load){
      ## Build RHugin domains
      res <- buildMixtureDomains(n.unknown, data, domainnamelist = domainnamelist,
                                 write = write, dir = dir, ntraces = ntraces,
                                 triangulate = triangulate, compile = compile, compress = compress, use.order = use.order)

      ## store the compression percentages
      out$compressionpct <- res$compressionpct

      ## list of pointers to domains
      out$domains <- res$domains
    }
    else {
        stop("Loading rhd-files is not available in DNAmixturesLite")
      ## read the rhd-files and store the list of pointers
      #out$domains <- lapply(domainnamelist, function(d){read.rhd(paste(dir, d, ".hkb", sep = ""))})
    }

    ## Add the names for binary nodes as attributes to each domain
    out$domains <- lapply(out$markers, function(m){
      als <- seq_len(nrow(out$data[[m]])) ## 1:(number of alleles)
      height <- out$heightnames
      O <- D <- Q <- list()
      for (r in 1:ntraces){
        if (r %in% observed){
          O[[r]] <- paste("O",r, als, sep = "_")
          D[[r]] <- paste("D",r, als, sep = "_")
          Q[[r]] <- paste("Q",r, als, sep = "_")
        } else {
          O[[r]] <- D[[r]] <- Q[[r]] <- character(0)
        }
      }
      O <- lapply(1:ntraces, function(r)if(all(is.na(out$data[[m]][,r+1])))NULL else paste("O",r, als, sep = "_"))
      D <- lapply(1:ntraces, function(r)if(all(is.na(out$data[[m]][,r+1])))NULL else paste("D",r, als, sep = "_"))
      Q <- lapply(1:ntraces, function(r)if(all(is.na(out$data[[m]][,r+1])))NULL else paste("Q",r, als, sep = "_"))

      structure(out$domains[[m]],
                O = O,
                D = D,
                Q = Q,
                n = outer(als, seq_len(n.unknown),
                  function(a, i)paste("n", i, a, sep = "_")))
    })
    names(out$domains) <- out$markers
  }

  else {
    ## no unknown contributors
    out$U <- character(0)
  }

  class(out) <- "DNAmixture"
  out
}



##' @rdname DNAmixture
##' @param x An object of class \code{DNAmixture}.
##' @param ... not used.
##' @export
##' @method print DNAmixture
print.DNAmixture <- function(x, ...){

  ## Printing the model:
  cat("\nA DNA mixture model with", x$k, "contributors.\n\n")
  if (length(x$K)){
    cat("Known:", x$K, "\n")
  }
  else cat("No known contributors\n")

  if (x$n.unknown > 0){
    cat("Unknown:", x$U, "\n\n")
  }
  else cat("No unknown contributors\n\n")

  cat("Mixtures included:", x$call, "\n")
  cat("Detection threshold(s):", unlist(x$C), "\n\n")
  invisible(x)
}
