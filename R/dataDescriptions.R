##' @name MC15
##' @title The MC15 DNA mixture
##' @description Peak heights from analysis of a bloodstain with DNA
##' from an unknown number of contributors. There are three reference
##' profiles \code{K1}, \code{K2}, and \code{K3} for individuals
##' related to the case.
##' @docType data
##' @format A \code{data.frame}.
##' @source Peter Gill et. al. (2008)
##' \emph{Interpretation of complex {DNA} profiles using empirical models and a method to measure their robustness}.
##' Forensic Science International: Genetics, 2(2):91--103.
NULL

##' @name MC18
##' @title The MC18 DNA mixture
##' @description Peak heights from analysis of a bloodstain with DNA
##' from an unknown number of contributors. There are three reference
##' profiles \code{K1}, \code{K2}, and \code{K3} for individuals
##' related to the case.
##' @docType data
##' @format A \code{data.frame}.
##' @source Peter Gill et. al. (2008)
##' \emph{Interpretation of complex {DNA} profiles using empirical models and a method to measure their robustness}.
##' Forensic Science International: Genetics, 2(2):91--103.
NULL

##' @name USCaucasian
##' @title The data base of allele frequencies for 302 US Caucasian profiles.
##' @docType data
##' @description Allele frequencies for US Caucasians.
##' @note There are a few errata found after publication of the allele
##' frequencies, but we have not updated the allele frequencies
##' accordingly. The frequencies also did not add up to 1, and this is
##' simply corrected by normalising within each marker.
##' @format A \code{data.frame}.
##' @source \url{http://www.cstl.nist.gov/strbase/NISTpop.htm}
NULL

##' @name UKCaucasian
##' @title Allele frequencies for UK Caucasians
##' @docType data
##' @description Database of allele frequencies for UK Caucasians.
##' @format A \code{data.frame}
##' @source The frequencies are produced from the allele counts for
##' one (the UK Caucasian) of three British sub-populations as found
##' on David Balding's homepage under his
##' \href{https://sites.google.com/site/baldingstatisticalgenetics/software/likeltd-r-forensic-dna-r-code}{\pkg{likeLTD}}
##' software for mixture analysis in \R.
##' @references \url{https://sites.google.com/site/baldingstatisticalgenetics/software/likeltd-r-forensic-dna-r-code}
NULL

##' @name NGM
##' @title NGM allele frequencies
##' @docType data
##' @description NGM allele frequencies. 
##' @format A list containing 
##' \describe{
##' \item{\code{USCaucasian}}{A \code{data.frame} with, for each marker,
##' the set of possible alleles as well as their corresponding frequency.}
##' }
##' @source Budowle et. al (2011) \emph{Population Genetic Analyses of the NGM STR loci}.
##' International Journal of Legal Medicine, 125(1):101-109.
NULL

##' @name SGMplusDyes
##' @title Dyes used for SGMplus
##' @docType data
##' @description Dyes used for the AmpFlSTR SGM Plus PCR Amplification Kit
##' @source Applied Biosystems
##' @format A list with one item per dye, each containing a vector of
##' marker-names specifying the order in which they occur in an EPG.
NULL

##' @name ProfilerDyes
##' @title Dyes used for Profiler plus
##' @docType data
##' @description Dyes used for the AmpFlSTR Profiler Plus PCR Amplification Kit
##' @source Applied Biosystems
##' @format A list with one item per dye, each containing a vector of
##' marker-names specifying the order in which they occur in an EPG.
NULL

##' @name NGMDyes
##' @title Dyes used for NGM
##' @docType data
##' @description Dyes used for the AmpFlSTR NGM PCR Amplification Kit
##' @source Applied Biosystems
##' @format A list with one item per dye, each containing a vector of
##' marker-names specifying the order in which they occur in an EPG.
NULL
