##' Create a data set for a DNA mixture model
##'
##' @description The function is intended for internal use in
##' \code{\link{DNAmixture}} and its purpose is to combine peak height
##' data, any reference profiles, and allele frequencies.
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
##' @param data Either one \code{data.frame} containing variables
##' \code{markers} and \code{allele}, \code{height}, or a list of one
##' or more of such \code{data.frame}s, corresponding to each EPG. If
##' a \code{data.frame} with reference profiles are not specified
##' separately, these should be found amongst the datasets.
##'
##' @param K Names for reference profiles.
##'
##' @param database A \code{data.frame} with variables \code{marker},
##' \code{allele}, and allele \code{frequency}.
##'
##' @param reference.profiles Optionally, a \code{data.frame} containing allele counts for
##' each reference profile specified in \code{K}. There should be one variable for each name
##' in \code{K} as well as variables \code{marker} and \code{allele}.
##'
##' @return list of \code{data.frame}'s indexed by markers and
##' containing variables
##' \item{\code{marker}}{The STR marker}
##' \item{\code{allele}}{Repeat number of the allele}
##' \item{\code{frequency}}{Allele frequency}
##' \item{\code{height1},\ldots, \code{heightN}}{Peak heights for each of the N mixtures analysed;
##' their order follows the order in which the mixtures are specified in \code{data}.}
##' \item{\code{stutter.from}}{For internal use. Where do the alleles get stutter from?
##' A value of \code{-1} means that the allele cannot receive stutter, corresponding to \code{gets_stutter} being false.
##' \code{allele[i]} receives stutter from \code{allele[stutter.from[i]]}}
##' \item{stutter.to}{As above. \code{allele[i]} can stutter to \code{allele[stutter.to[i]]}}
##' as well as known DNA profiles as labelled by \code{K}.
##' @examples
##' ## Create a dataset for two markers with each 3 observed alleles
##' epgdf <- data.frame(marker = rep(c("FGA", "TH01"), each = 3),
##'                     allele = c(18, 23, 27, 7, 8, 9.3), ## Observed alleles
##'                     height = c(100, 100, 200, 200, 100, 100), ## Peak heights
##'                     Anna = c(0,0,2,1,0,1),             ## Anna's profile
##'                     Peter = c(1,1,0,1,1,0))            ## Peter's profile
##'
##' data(USCaucasian)
##' dat <- DNAmixtureData(epgdf, K = c("Anna", "Peter"), database = USCaucasian)
##' @export
##' @author Therese Graversen
DNAmixtureData <- function(data, database, K = character(0), reference.profiles = NULL) {

  if(is(data, "data.frame"))
    dflist <- list(data)
  else
    dflist <- data

  ## All markers considered in the collection of traces.
  ## (so extra markers in reference profiles or database will be dropped)
  allmarkers <- unique(unlist(lapply(dflist, "[[", "marker")))
  K.dat <- if(is.null(reference.profiles)) K else character(0)

  dflist <- lapply(dflist, function(d){
    if(!(class(d$allele) %in% c("numeric", "integer"))){
      stop("Alleles should be numerical values")
    }
    out <- d[,names(d) %in% c("marker", "allele", "height", K.dat)]
    ## Integer code for alleles
    out$allele <- as.integer(round(d$allele*10, digits = 0))
    out
    })

  for (r in seq_along(dflist)){
    ## Rename height variable for each dataset; index corresponds to EPG.
    names(dflist[[r]])[names(dflist[[r]])=="height"] <- paste0("height",r)
  }

  dat <- Reduce(function(x, y){
    known.names <- intersect(K.dat, intersect(names(x), names(y)))
    merge(x, y, by = c("marker", "allele", known.names), all = TRUE)},
                dflist)

  if (!is.null(reference.profiles)){
    if (!all(K %in% names(reference.profiles)))
      stop ("Reference profiles do not contain the specified contributors")
    if(!(class(reference.profiles$allele) %in% c("numeric", "integer")))
      stop("Allele names in the reference profiles are not numerical")

    ref <- reference.profiles[reference.profiles$marker %in% allmarkers, , drop=FALSE]
    ref <- ref[, names(ref) %in% c("marker", "allele", K), drop=FALSE]
    
    ref$allele <- as.integer(round(ref$allele*10, digits = 0))
    dat <- merge(dat, ref, by = c("marker", "allele"), all = TRUE)
  }

  ## merge on the allele frequencies.
  if(!(class(database$allele) %in% c("numeric", "integer"))){
    stop("Allele names in the database are not numerical")
  }
  db <- database[database$marker %in% allmarkers,]
  db$allele <- as.integer(round(db$allele*10, digits = 0))
  dat <- merge(dat, db, by = c("marker", "allele"), all=TRUE)

  if(any(is.na(dat$frequency))){
    warning("Some alleles in the data are not found in the database and are assigned frequency NA", call. = FALSE, immediate. = TRUE)
  }

  ## Sort variables as allele, height1 ,..., heightR, K, marker, and the allele integer code
  vars <- c("allele", paste0("height",seq_along(dflist)), K, "frequency", "marker")
  ## othervars <- names(dat)[!(names(dat) %in% vars)]

  dat <- dat[, vars, drop = FALSE]

  ## Split into a data.frame per marker
  ## (and drop any non-used levels of marker)
  dat <- split(dat, dat$marker, drop = TRUE)

  ## Now deal with any NA in heights or reference profiles.
  prepare.marker <- function(d){

    ## If no heights observed for a trace at this marker, leave NA's.
    ## If some heights are observed, set remaining NA's to 0.
    ## *problem here if kits differ in allele range*
    for (v in seq_along(dflist)+1){
      if (!all(is.na(d[,v]))){
        d[,v] <- ifelse(is.na(d[,v]), 0, d[,v])
      }
    }

    ## Check that reference profiles are specified for this marker.
    for (v in K){
      if (all(is.na(d[,v])))
        stop("Reference profiles are not available for all markers")
      else d[,v] <- ifelse(is.na(d[,v]), 0, d[,v])
      if (sum(d[,v]) != 2){
        stop("All contributors should have exactly 2 alleles per marker")
      }
    }

    ##    if(any(colSums(subset(d, select = K)) != 2))
    ##      stop("All contributors should have exactly 2 alleles per marker")

    ## We do not want markers as factors (or at least then we should drop unused factor levels)
    d$marker <- as.character(d$marker)

    ## Make sure alleles are numeric, and add on stutter (except for amelogenin)
    if (unique(d$marker) == "AMEL"){
      ## ordering alleles as 0(X), 1(Y); removing rownames.
      d <- d[order(d$allele),]
      row.names(d) <- NULL

      ## No stutter for Amelogenin
      d$stutter.from <- -1
      d$stutter.to <- -1
      d$gets_stutter <- FALSE
      d$can_stutter <- FALSE
    }
    else {
      ## ordering alleles according to minus one stutter
      ## d <- d[order(d$allele - floor(d$allele), d$allele),]
      d <- d[order(d$allele %% 10, d$allele),]

      ## Remove the rownames from the old ordering
      row.names(d) <- NULL

      ## Index of allele a' stuttering to allele a
      ## If allele[a]+1 exists in the set of alleles, then allele[a] can receive stutter
      d$stutter.from <- match(d$allele + 10L, d$allele, nomatch = -1)

      ## Index of allele a' receiving stutter from allele a
      ## If allele[a]-1 exists in the set of alleles, then allele[a] can stutter
      d$stutter.to <- match(d$allele - 10L, d$allele, nomatch = -1)

      ## Booleans indicating whether an allele can receive stutter or stutter itself.
      d$gets_stutter <- d$stutter.from != -1
      d$can_stutter <- d$stutter.to != -1
    }
    ## Allele back to decimal format
    d$allele <- d$allele/10
    d
  }

  lapply(dat, prepare.marker)
}
