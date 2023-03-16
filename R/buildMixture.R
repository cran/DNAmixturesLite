##' Create RHugin domains for computation.
##'
##' @description The function is intended for internal use in
##' \code{\link{DNAmixture}}, and it creates one network per
##' marker.
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
##' For one marker, the network contains variables
##' \describe{
##' \item{\code{n_i_a}}{The number of alleles \eqn{a} that unknown contributor \code{Ui} possesses.}
##' \item{\code{S_i_a}}{Cumulative allele counts, summarising how many alleles amongst types 1, \ldots, \eqn{a} that contributor \eqn{i} possesses.}
##' \item{\code{O_r_a}}{Binary variable representing observed peak height for allele \eqn{a} in EPG \eqn{r}. See \code{\link{setCPT.O}} for details.}
##' \item{\code{D_r_a}}{Binary variable representing the event of a peak for allele \eqn{a} falling below threshold in EPG \eqn{r}.}
##' \item{\code{Q_r_a}}{Binary variable representing the event of a peak for allele \eqn{a} being smaller than the peak observed in EPG \eqn{r}.}
##' }
##'
##' The network is by default triangulated, compiled,
##' compressed, and, optionally, saved as a hkb-file. If the networks
##' are large -- if there are many unknown contributors -- it is worth
##' considering saving the networks rather than re-building them every time
##' \code{DNAmixture} is called.
##'
##' @note If amelogenin is included as a marker, use
##' integer codes for the alleles with X preceding Y, e.g. X=0, Y=1.
##'
##' @param n.unknown Number of unknown contributors
##' @param data A list of \code{data.frame}s on the form returned by
##' \code{\link{DNAmixtureData}}.
##' @param domainnamelist List of names for networkfiles.
##' @param write [not available in lite-version] Save networkfiles? Defaults to \code{FALSE}.
##' @param dir Path to the networkfiles.
##' @param ntraces Number of traces (EPG's)
##' @param triangulate Triangulate the networks? Default is to
##' triangulate the network using a good elimination order.
##' @param compile Compile the networks? Defaults to \code{TRUE}.
##' @param compress Compress the network? Defaults to \code{TRUE} and
##' is strongly recommended for models with a large number of
##' contributors.
##' @param use.order Should the default elimination order be used for
##' triangulation?  Otherwise the "total.weight" heuristic for
##' triangulation in HUGIN is used.
##'
##' @return A list containing a list of pointers for RHugin
##' domains. Additionally a list of the total size of the compressed
##' junction tree in percent of the uncompressed size.
##' @author Therese Graversen
##' @export
buildMixtureDomains <- function(n.unknown, data, domainnamelist, write, dir = NULL, ntraces = 1,
                                triangulate = TRUE, compile = TRUE, compress = TRUE, use.order = TRUE){

  ## function for creating network for one marker
  .buildMixtureDomain <- function(domainname, n.unknown, data, write, dir, AMEL = FALSE){

    ## Build the network:
    domain <- hugin.domain()
    compressionpct <- NULL

    ## Extract number of alleles, Frequencies of alleles and indices for
    ## stutter contributions.
    alleles <- seq_along(data$allele)
    allele.freqs <- data$frequency
    stutter.from <- data$stutter.from

    ## Matrices of names for nodes: n_ia, S_ia and binary nodes O_ra etc.
    n <- outer(1:n.unknown, alleles, function(i,a)paste("n", i, a, sep = "_"))
    S <- outer(1:n.unknown, alleles, function(i,a)paste("S", i, a, sep = "_"))
    O <- outer(1:ntraces, alleles, function(r,a)paste("O", r,a, sep = "_"))
    D <- outer(1:ntraces, alleles, function(r,a)paste("D", r,a, sep = "_"))
    Q <- outer(1:ntraces, alleles, function(r,a)paste("Q", r,a, sep = "_"))

    ## local function for creating all n_ia and S_ia for one contributor
    add.contributor <- function(i){
      add.allele.for.contributor <- function(a){

        ## Create a node corresponding to n_ia
        n_ia <- n[i, a]
        add.node(domain, n_ia, states = 0:2, subtype = "numbered")

        ## Create nsum_iAllele(a) = nsum_i{Allele(a-1)} + n_i{Allele(a-1)}
        S_ia <- S[i, a]
        add.node(domain, S_ia, states = 0:2, subtype = "numbered")
        add.edge(domain, S_ia, n_ia)

        if (a == 1){
          ## n_ia is binomial(2, allele.freqs[a])
          tab <- get.table(domain, n_ia)
          ns <- tab[,n_ia]
          if (AMEL){
            ## n_X ~ 1 + bin(1, q_X)
            tab$Freq <- dbinom(ns - 1, 1, allele.freqs[a]/sum(allele.freqs))
          }
          else {
            ## All other markers
            tab$Freq <- dbinom(ns, 2, allele.freqs[a]/sum(allele.freqs))
          }
          set.table(domain, n_ia, tab)

          ## nsum_ia = n_ia
          tab <- get.table(domain, S_ia)
          ns <- tab[,n_ia]
          nsums <- tab[,S_ia]
          tab$Freq <- ifelse(nsums == ns, 1, 0)
          set.table(domain, S_ia, tab)
        }

        if (a > 1) {
          ## n_ia|nsum_i(a-1) is binomial(2 - nsum_i(a-1), sum_{j=a}^A q_j)
          ## provided that the sum of alleles is at most 2.
          ## For a sum > 2 we set n_ia = 0.
          ## Note, this only uses the ordering of alleles, not the stutter-ordering

          S_i_prev <- S[i, a-1]
          add.edge(domain, n_ia, S_i_prev)

          tab <- get.table(domain, n_ia)
          prev.nsums <- tab[, S_i_prev]
          ns <- tab[, n_ia]
          s <- sum(tail(allele.freqs, -(a-1))) #sum(allele.freqs[a:length(allele.freqs)])
          ## Need to make sure that we do not divide by 0,
          ## when the remaining allele frequencies are all 0.
          tab$Freq <- dbinom(ns, 2 - prev.nsums,
                             allele.freqs[a]/ifelse(s>0,s,1)
                             )
          set.table(domain, n_ia, tab)

          ## nsum_ia = nsum_i(a-1) + n_ia
          add.edge(domain, S_ia, S_i_prev)

          tab <- get.table(domain, S_ia)
          nsums <- tab[, S_ia]
          prev.nsums <- tab[, S_i_prev]
          ns <- tab[, n_ia]
          allelesums <- ns + prev.nsums

          ## For valid combinations of parents (sum to 2 or less), the node is a
          ## sum of the parents. For other combinations (have probability 0) we
          ## use a uniform distribution
          tab$Freq <- ifelse(allelesums <= 2,
                             ifelse(nsums == allelesums, 1, 0),
                             1/3)

          set.table(domain, S_ia, tab)
        }
      }
      lapply(alleles, add.allele.for.contributor)
    }

    ## local function for creating O_a, D_a and Q_a
    addBinaryNodes <- function(a){

      ## parent set is n_ia and n_ia', if an a' can stutter to a
      ind <- c(a, stutter.from[[a]][stutter.from[[a]] != -1])
      parents <- as.vector(n[1:n.unknown, ind])

      for (r in 1:ntraces){
        ## Binary nodes
        add.node(domain, O[r,a], subtype = "boolean") ## for likelihood function
        add.node(domain, D[r,a], subtype = "boolean") ## for P(height < C)
        add.node(domain, Q[r,a], subtype = "boolean") ## for P(height < obs_height)

        ## Add an edge from each parent to the binary node
        lapply(parents, function(p)add.edge(domain, O[r,a], p))
        lapply(parents, function(p)add.edge(domain, D[r,a], p))
        lapply(parents, function(p)add.edge(domain, Q[r,a], p))

        ## Note: CPT's for binary nodes are set to the default uniform distributions
        ## when creating the network.
      }
    }

    ## Create all n_ia and S_ia
    lapply(1:n.unknown, add.contributor)

    ## Adding the binary nodes
    lapply(alleles, addBinaryNodes)

    ## Triangulate using either the total.weight heuristic, or a specified triangulation
    ## The default triangulation is optimal when all alleles are in one "stutter-group"
    if (triangulate) {
      if (use.order){
        A <- length(alleles)
        if (AMEL)
          o <- c(O, D, Q, S[,c(A, A-1)], n)
        else
          o <- c(O, D, Q, S[,c(A, A-1, 1)], n[,1], rbind(n[,2:(A-2)], S[,2:(A-2)]), n[,c(A-1,A)])
        triangulate(domain, order = o)
      }
      else {
        triangulate(domain, method = "total.weight")
      }
      ## compile
      if (compile) {
        compile(domain)
        ## compress and save the relative sizes
        if (compress) compressionpct <- compress(domain)
      }
    }

    ## Write network to file, if desired. -- not available for lite-version
      if (write){
          warning("Writing rhd files is not available in DNAmixturesLite")
      ## Create directories if they do not exist
      #dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      #write.rhd(domain, filename = paste(dir, domainname, ".hkb", sep = ""), type = "hkb")
    }

    list(compressionpct = compressionpct,
         domain = domain)
  }

  ## Iterate over all markers
  markers <- names(data)

  li <- lapply(markers, function(m){
    .buildMixtureDomain(domainname = domainnamelist[[m]], data = data[[m]],
                        n.unknown = n.unknown, write = write, dir = dir, AMEL = (m == "AMEL"))
  })
  names(li) <- markers

  list(
    ## Percent of original size
    compressionpct = lapply(li, function(i){i$compressionpct}),
    ## domain pointers
    domains = lapply(li, function(i){i$domain})
    )
}
