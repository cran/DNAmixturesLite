##' @importFrom grDevices dev.interactive devAskNewPage
##' @importFrom graphics abline axis box clip curve lines mtext par plot.new points
##' @importFrom methods is
##' @importFrom stats dbinom dgamma pgamma ppoints predict qnorm rgamma simulate
##' @importFrom utils data head tail
##' @importFrom Rsolnp solnp
##' @importFrom Matrix bdiag
##' @importFrom numDeriv hessian
##' @importFrom gRaven add.edge add.node get.belief get.normalization.constant get.table hugin.domain initialize.domain map.configurations retract set.finding set.table propagate.gRaven
##' @importFrom gRbase propagate
NULL


##' Statistical Inference for Mixed Samples of DNA (Lite-Version)
##'
##' @description Tools for statistical inference for one or multiple DNA mixtures.
##'
##' \emph{IMPORTANT: This is the \pkg{DNAmixturesLite} package, which is intended as a service to enable users to try \pkg{DNAmixtures} without purchasing a commercial licence for Hugin. When at all possible, we strongly recommend the use of \pkg{DNAmixtures} rather than this lite-version. See \url{https://dnamixtures.r-forge.r-project.org/} for details on both packages.}
##'
##'  \emph{While the lite-version seeks to provide the full functionality of \pkg{DNAmixtures}, note that computations are much less efficient and that there are some differences in available functionality. Be aware that the present documentation is copied from \pkg{DNAmixtures} and thus may not accurately describe the implementation of this lite-version.}
##'
##' @details The package implements a statistical model for analysis
##' of one or more mixed samples of DNA in the possible presence of
##' dropout and stutter. Details of the model can be found in Cowell
##' et. al (2013), and details on the model checking tools and
##' Bayesian network structure can be found in Graversen and Lauritzen
##' (2014).
##'
##' Any hypothesis involving unknown contributors relies on
##' computations in a Bayesian network. For performing such
##' computations, \pkg{DNAmixtures} package relies on Hugin
##' (\url{https://www.hugin.com}) through the \R-package \pkg{RHugin}. For an
##' installation guide, see the package webpage
##' \url{https://dnamixtures.r-forge.r-project.org}.
##'
##' Although \pkg{DNAmixtures} can be installed with only the free
##' version of Hugin, the size of the networks will in practice
##' require the full licence. In theory, the implementation allows
##' analysis with an arbitrary number of unknown
##' contributors. However, in practice, depending on hardware and
##' time-constraints working with up to 5 or 6 unknown contributors
##' seems realistic.
##'
##' @section Summary of the statistical model:
##'
##' The statistical model jointly models the observed peak heights and
##' the set of contributors to the DNA sample(s). In the event of
##' analysing multiple DNA mixtures, the union of the contributors is
##' used as the contributor set for each mixture. By allowing a
##' contribution of zero, we cover the case of a contributor not
##' having contributed to a particular mixture.
##'
##' Genotypes for unknown contributors are modelled using
##' allele-frequencies from a database specified by the user. The
##' database is also used to define the range of alleles at each
##' marker. A genotype for an unknown contributor is represented by a
##' vector of allele counts \eqn{n_{ia}}, counting for each allele
##' \eqn{a} the number of alleles \eqn{i} that a person possesses; in
##' the network for a marker, the allele count \eqn{n_{ia}} is
##' represented by a variable \code{n_i_a}. The vector of allele
##' counts follows a multinomial distribution with \eqn{\sum_i n_{ia}
##' = 2} and the specified allele frequencies. It is assumed that
##' genotypes are independent across markers and between
##' contributors. If desired, the database of allele frequencies may
##' be corrected for F_st or sampling adjustment before use.
##'
##' Peak heights are assumed mutually independent and their
##' distributions for a fixed set of DNA profiles are modelled using
##' gamma distributions.  The peak height for allele \eqn{a} in EPG \eqn{r}
##' is assumed to follow a gamma distribution with scale parameter
##' \eqn{\eta_r} and shape parameter
##'
##' \deqn{\rho_r \sum_a ((1-\xi_{ra})n_{ia} + \xi_{r,a+1} n_{i,a+1})\phi_{ri}.}
##'
##' Applying a detection threshold \eqn{C_r\ge 0}, any peak height
##' falling below the threshold is considered to be 0.  The peak
##' heights are denoted by \code{height1, \ldots, heightR}.
##'
##' The model parameters are for each DNA mixture
##' \describe{
##' \item{\eqn{\phi}}{The proportions of DNA from each contributor.}
##' \item{\eqn{\rho}}{Amplification parameter, which will be larger for larger amounts of DNA amplified.}
##' \item{\eqn{\eta}}{Scale parameter for the gamma distribution.}
##' \item{\eqn{\xi}}{Mean stutter percentage. Allele \eqn{a} uses stutter parameter \eqn{\xi_a = \xi} if the allele \eqn{a-1} is included in the model, and \eqn{\xi_a = 0} otherwise}
##' }
##'
##' An alternative parametrisation uses \eqn{\mu = \rho \eta} and
##' \eqn{\sigma = 1/\sqrt{\rho}}, which can be interpreted as the mean
##' peak height and the coefficient of variation respectively. Besides
##' being interpretable, an advantage of this reparametrisation is
##' that the parameters are fairly orthogonal.
##'
##' The model assumes the model parameters to be the same across
##' markers.  Relaxations of these assumptions are not implemented
##' here.
##'
##'
##' @section Computation by auxiliary variables:
##'
##' The computational approach of the implementation of this package
##' is discussed in Graversen and Lauritzen (2014).
##'
##' The Bayesian networks include three types of auxiliary variables
##' \code{O}, \code{D}, and \code{Q}; these can be thought of as
##' representing the observed peak heights, the absence/presence of
##' peaks, and the peak height distribution function. Note that if
##' invalid tables are set -- for instance if very extreme parameter
##' values are used, or if the vector of mixture proportions is
##' mis-labeled -- then any subsequent propagation will fail. No
##' roll-back functionality has so far been implemented to fix this,
##' and the easiest solution is to re-fit the mixture model.
##'
##' The workhorses of this package are the functions
##' \code{\link{setCPT.O}}, \code{\link{setCPT.D}} and
##' \code{\link{setCPT.Q}} for setting the conditional probability
##' tables for the three types of auxiliary variables according to
##' specified peak heights and model parameters.
##'
##' @section Amelogenin: As an experiment, it is possible to add the
##' marker Amelogenin, provided that the marker is named "AMEL" and
##' that the coding of alleles X and Y is of a particular form. One
##' example of a suitable form is the coding X = 0 and Y = 1. The
##' allele frequencies used should then also contain a marker "AMEL",
##' and here frequencies have a slightly different interpretation than
##' for the rest of the markers; as all people possess one X, the
##' frequencies of X and Y denote the presence of an additional X or Y
##' respectively, and thus the frequencies correspond directly to the
##' proportions of the two genders.
##'
##' @name DNAmixturesLite-package
##' @aliases DNAmixturesLite DNAmixtures
##' @docType package
##'
##' @author Therese Graversen \email{theg@@itu.dk}
##'
##' @references Details on the implemented model may be found in
##'
##' Cowell, R. G., Graversen, T., Lauritzen, S., and Mortera, J. (2015).
##' \emph{Analysis of Forensic DNA Mixtures with Artefacts}. With supplementary material documenting the analyses using \pkg{DNAmixtures}.
##' Journal of the Royal Statistical Society: Series C (Applied Statistics).
##' Volume 64, Issue 1, pages 1-48.
##'
##' Graversen, T. (2014)
##' \emph{Statistical and Computational Methodology for the Analysis of Forensic DNA Mixtures with Artefacts}.
##' DPhil. University of Oxford.
##' \url{https://ora.ox.ac.uk/objects/uuid:4c3bfc88-25e7-4c5b-968f-10a35f5b82b0}.
##'
##' Graversen, T. and Lauritzen, S. (2014).
##' \emph{Computational aspects of {DNA} mixture analysis}.
##' Statistics and Computing, DOI: 10.1007/s11222-014-9451-7.
##'
##' @example inst/examples/main.R
##'
NULL
