##' Plot a DNA mixture model
##'
##' @description Plot of peak heights against repeat numbers for each marker.
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
##' @export
##' @method plot DNAmixture
##' @param x A \code{DNAmixture}.
##' @param traces Indices giving the mixtures, for which plots should be made.
##' @param markers Vector of names giving the markers, for which plots should be made.
##' @param epg Should a stylized EPG be produced? This requires a list
##' of dyes to be specified for the mixtures, possibly through
##' \code{\link{dyes}}.
##' @param dyecol List containing for each mixture a vector of dye
##' names. The default is to use the names in \code{dyes(x)}. Set to
##' \code{NULL} to ignore.
##' @param pw Peaks are \code{2*pw} wide.
##' @param threshold Should the detection thresholds be indicated on the plot?
##' @param add Add to existing plot? (not meaningful if \code{EPG = TRUE}.)
##' @param ask Should the user be prompted for page changes? The
##' default is \code{TRUE}, when the device is
##' \code{dev.interactive()} and there is a need for multiple pages.
##' @param panel Alternative function for drawing the peaks. For
##' instance \code{lines} can be used for making triangular peaks. This functionality will be removed.
##' @param ... Other parameters to be passed on to \code{plot}
##' @return A \code{data.frame} containing the plotted data points.
##' @examples
##' data(MC15, MC18, USCaucasian, SGMplusDyes)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50,50), k = 3, K = c("K1", "K3"),
##'                   database = USCaucasian)
##' ## Plot as a stylized EPG, using "orange" for the "yellow" dye
##' names(SGMplusDyes) <- c("blue", "green", "orange")
##' dyes(mix) <- list(SGMplusDyes, SGMplusDyes)
##' plot(mix, epg = TRUE)
##' ## We can also supress the dye colors
##' plot(mix, traces = 1, epg = TRUE, dyecol = NULL)
##' ## Create a user specified layout
##' op <- par(mfrow = c(2,3))
##' plot(mix, markers = c("VWA", "D19S433", "TH01"), col = "red")
##' par(mfrow = op)
##' plot(mix, traces = 1, markers = "D19S433", col = "orange",
##'      main = "A single marker", cex.main = 2, ylim = c(0,200))
plot.DNAmixture <- function(x, traces = seq_len(x$ntraces), markers = x$markers,
                            epg = FALSE, dyecol = lapply(dyes(x), names),
                            pw = 0.4, threshold = TRUE, panel = NULL,
                            add = FALSE, ask = NULL, ...){

    if (is.null(ask)){
        ask <- dev.interactive() &&
            ( (!epg &&(length(markers)*length(traces) > prod(par("mfrow"))))|| (epg && length(traces)>1) )
    }

    ## Mixtures analysed
    if (any(traces > x$ntraces))stop("The specified traces are invalid")

    ## Largest peak amongst those we are plotting
    maxheight <- max(sapply(x$data[markers], function(d)max(d[,traces+1, drop = FALSE])), na.rm = TRUE)

    ## Setting and restoring the interaction with the user
    if (ask) {
        oldask <- devAskNewPage(ask)
        on.exit(devAskNewPage(oldask))
    }

    if(epg){
        ## Need to change, and therefore restore graphical parameters
        oldpar <- par(oma = c(0,0,3,0),         # No outer margin (inches)
                      mar = c(2.5,2.5,1.5,0.5), # margin between plot- and figure-region. (lines)
                      tcl = -0.2,               # tick-length (1=lineheight)
                      mgp = c(1.3,.3,0))        # position of axislabel, ticklabels and axis (lines)
        on.exit(par(oldpar), add = TRUE)
    }

    ## Main plot function to be called at the end
    plotmix <- function(x, xlim = NULL, ylim = NULL, xlab = "Allele", ylab = "Peak height",
                        col = par("col"), bg = NA, pch = par("pch"),
                        cex = par("cex"), lty = par("lty"), lwd = par("lwd"),
                        axes = TRUE, frame.plot = axes, xaxt = par("xaxt"), yaxt = par("xaxt"),
                        ann = par("ann"), cex.lab = par("cex.lab"),
                        col.lab = par("col.lab"), font.lab = par("font.lab"),
                        cex.axis = par("cex.axis"), col.axis = par("col.axis"),
                        font.axis = par("font.axis"),
                        main = NULL, col.main = par("col.main"), font.main = par("font.main"), cex.main = par("cex.main"),
                        ...){

        ## Function for peak heights and positions
        peakPoints <- function(alleles, heights){
            ## Treating all markers as having 4 bp repeats.
            alleles <- floor(alleles) + 2.5* (alleles - floor(alleles))
            xs <- as.vector(outer(c(-pw/4,0,pw/4), alleles, "+"))
            ys <- as.vector(rbind(0,heights,0))
            o <- order(xs)
            list(x = xs[o], y = ys[o], mid = alleles, height = heights)
        }
        ## The plot information for *all mixtures*, not just the ones we plot
        pp <- lapply(x$data[markers],
                     function(d){apply(d[,seq_len(x$ntraces)+1, drop = FALSE], 2,
                                       function(height)peakPoints(d$allele, height))})

        ## Curves for peak plotting
        nicePeak <- function(a,b, ...){
            f <- function(t)3*t^2-2*t^3
            for(i in 1:(length(a)-1)){
                x1 <- a[i:(i+1)]
                y1 <- b[i:(i+1)]
                curve(y1[1]*(1-f((x-x1[1])/(x1[2]-x1[1]))) + y1[2]*f((x-x1[1])/(x1[2]-x1[1])), x1[1],x1[2], add = TRUE, ...)
            }
        }

        ## What to plot for a combination of marker m and mixture r
        plotFun <- function(m, r, col = col, ...){

            if (!add){ ## Start a new plot
                ylim <- if (is.null(ylim)) c(0,maxheight) else ylim
                xlim <- if (is.null(xlim)) range(pp[[m]][[r]]$x) else xlim
                main <- if (is.null(main)) m else main
                plot(0,0, type = "n",
                     xlab = xlab, ylab = ylab,
                     main = main, ylim = ylim, xlim = xlim, axes = FALSE,
                     col.main = col.main, font.main = font.main, cex.main = cex.main)

                if (axes){
                    ## Draw axes
                    a <- x$data[[m]]$allele
                    scaled.a <- pp[[m]][[r]]$mid
                    axis(1, at = scaled.a[a-floor(a) == 0], labels = as.integer(a[a-floor(a) == 0]), xaxt = xaxt, ...)
                    axis(1, at = scaled.a[a-floor(a) != 0], labels = FALSE, tcl = par("tcl")*.5, xaxt = xaxt, ...)
                    axis(2, yaxt = yaxt, ...)
                }
                if (frame.plot){
                    ## Draw box around the plot
                    box(...)
                }
            }

            ## Peaks are drawn by the panel function (default smooth peaks)
            panel <- if (is.null(panel)) nicePeak ## else match.fun(panel)
            do.call(panel, list(pp[[m]][[r]]$x, pp[[m]][[r]]$y, col = col, lwd = lwd, lty = lty, ...))

            ## Could annotate with known profiles --- not implemented

            ## Add the threshold line
            if (threshold){
                abline(h = x$C[[r]], lty = "22")
            }
        }

        ## Assembling everything
        if (epg){
            ## Need to change plot layout, so need to restore that on exit
            old.mfrow <- par("mfrow")
            on.exit(par(mfrow = old.mfrow), add = TRUE)
            for (r in traces){
                submarkers <- intersect(names(x$observed)[sapply(x$observed, function(o)r %in% o)], markers)

                if (is.null(x$dyes)){stop("Dyes need to be set to plot the EPG")}
                ## One row in the plot per dye
                dyes <- x$dyes[[r]]
                submarkers <- intersect(unlist(dyes), submarkers) ## should be in the order of dyes
                rs <- length(dyes) ## number of dyes
                cs <- max(sapply(dyes, length)) ## number of markers for a dye
                par(mfrow = c(rs, cs))
                for (i in 1:rs){
                    for (j in 1:cs){
                        if (j <= length(dyes[[i]]) && dyes[[i]][j] %in% submarkers){
                            plotFun(dyes[[i]][j], r, col = if (!is.null(dyecol))dyecol[[r]][i] else col,...)
                        } else plot.new()
                    }
                }
                mtext(paste("Mixture", r), line = 1, cex = cex.main, outer = TRUE)
            }
        }
        else{
            for(r in traces){
                ## Plot in the order of the specified markers
                submarkers <- intersect(markers, names(x$observed)[sapply(x$observed, function(o)r %in% o)])
                for(m in submarkers){plotFun(m, r, col = col, ...)}
            }
        }

        invisible(pp)
    }

    ## Call the main plot function
    plotmix(x, ...)
}

##' Quantile-Quantile plot for assessing the distribution of observed peak heights.
##'
##' @description
##' Given \eqn{Z_a \ge C}, the peak height \eqn{Z_a} follows a
##' continuous distribution. The probability transform \eqn{P(Z_a \le z_a|
##' Z_a \ge C)} thus follows a uniform distribution. The function
##' \code{qqpeak} computes the quantiles of the probability transform
##' and produces a quantile-quantile plot. Information about the other
##' peak heights in the EPG may be taken into account in the
##' distribution of a peak.
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
##' @export
##' @param x A \code{DNAmixture} model.
##' @param pars A \code{mixpar} parameter.
##' @param dist \code{"joint"}, \code{"conditional"}, or \code{"prequential"}. See \code{\link{predict.DNAmixture}} for details.
##' @param xlab Legend for the x-axis.
##' @param ylab Legend for y-axis.
##' @param xlim Range of x-axis. Default is a plot with equal ranges
##' on the x- and y-axis.
##' @param ylim Range of y-axis.
##' @param by.allele Order of conditioning when \code{dist="prequential"}. See \code{\link{predict.DNAmixture}}
##' @param plot Should a plot be produced? Defaults to \code{TRUE}.
##' @param ... Other arguments to be passed on to \code{plot}.
##' @return A \code{data.frame} containg quantiles \code{q}
##' corresponding to \eqn{P(Z < z_{obs} | Z > 0)} in the specified
##' distribution, together with other quantities computed by
##' \code{\link{predict.DNAmixture}}. The data are ordered according
##' to \code{q}.
##' @author Therese Graversen
##' @examples
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = c("K1", "K3"),
##'                   database = USCaucasian)
##' p <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08,0.08),
##'             phi = list(c(U1 = 0.1, K3 = 0.2, K1 = 0.7)))
##' qqpeak(mix, pars = p, dist = "conditional")
##' \donttest{
##' ## If desired, we can make a customised plot -- here we color according to the two mixtures
##' qq <- qqpeak(mix, pars = p, dist = "conditional", plot = FALSE)
##' plot(ppoints(qq$q), qq$q, xlab = "Uniform quantiles", ylab = "Probability transform",
##'      col = qq$trace, pch = qq$trace)
##' }
qqpeak <- function(x, pars, dist = c("joint", "conditional", "prequential"),
                   by.allele = TRUE, plot = TRUE,
                   xlab = "Uniform quantiles", ylab = NULL, xlim = NULL, ylim = xlim, ...){
    dist <- match.arg(dist)

    if (is.null(ylab)){
    ylab <- switch(dist,
                   joint = expression(paste("P(", Z[a] <= z[a], " | ", Z[a] >= C,")")),
                   conditional = expression(paste("P(", Z[a] <= z[a], " | ", Z[b] == z[b], ", ", b!=a, ", ", Z[a] >= C,")")),
                   prequential = expression(paste("P(", Z[a] <= z[a], " | ", Z[b] == z[b], ", ", b<a, ", ", Z[a] >= C,")"))
                   )
  }

  d <- do.call("rbind", predict(x, pars = pars, dist = dist, by.allele = by.allele))
  rownames(d) <- NULL
  d <- subset(d, d$height > 0)
  d$q <- (d$smaller - d$unseen)/d$seen

  d <- d[order(d$q), ]

  if (plot){
    xs <- ppoints(d$q)
    ys <- sort(d$q)

    if (is.null(xlim)){xlim <- range(xs, ys)}
    plot(xs, ys,
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, ...)
  }
  invisible(d)
}


##' Plot simulated peak heights for multiple DNA mixtures.
##'
##' @description A plot will be made for each combination of samples and markers
##' specified, and it is up to the user to specify a layout for the
##' plots (e.g. via calls to \code{\link{par}})
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
##' @param x A \code{DNAmixture}.
##' @param sims A set of simulated peak heights. See also \code{\link{rPeakHeight}}.
##' @param traces Selected traces to plot.
##' @param markers Selected markers to plot.
##' @param pw Peaks are \code{2*pw} wide.
##' @param ylim Range of the y-axis.
##' @param border Color of the boxes.
##' @param ... Arguments passed on to \code{plot.DNAmixture}.
##' @return Invisibly the peak heights as used for the boxplots.
##' @author Therese Graversen
##' @method boxplot DNAmixture
##' @export
##' @importFrom graphics boxplot
##' @examples
##'
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = c("K1", "K3"),
##'                   database = USCaucasian)
##' p <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = c(0.08,0.08),
##'             phi = list(c(U1 = 0.1, K3 = 0.2, K1 = 0.7), c(K1 = 0.9, K3 = 0.05, U1 = 0.05)))
##' rpm <- rPeakHeight(mix, p, nsim = 1000, dist = "joint")
##' oldpar <- par("mfrow")
##' par(mfrow = c(1,2))
##' boxplot(mix, rpm, traces = 1:2, markers = "VWA")
##' par(mfrow = oldpar)
##'
##' \donttest{
##' data(MC15, MC18, USCaucasian)
##' mix <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = "K3", database = USCaucasian)
##' p <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = c(0.08,0.08),
##'             phi = list(c(U2 = 0.1, K3 = 0.2, U1 = 0.7), c(U1 = 0.9, K3 = 0.05, U2 = 0.05)))
##' rpm <- rPeakHeight(mix, p, nsim = 50, dist = "joint")
##' rpc <- rPeakHeight(mix, p, nsim = 50, dist = "conditional")
##' oldpar <- par("mfrow")
##' par(mfrow = c(2,2))
##' boxplot(mix, rpm, traces = 1:2, markers = "VWA") ## First row of plots
##' boxplot(mix, rpc, traces = 1:2, markers = "VWA") ## Second row of plots
##' par(mfrow = oldpar)
##' mixK <- DNAmixture(list(MC15, MC18), C = list(50, 38), k = 3, K = c("K1", "K2", "K3"),
##'                     database = USCaucasian)
##' pK <- mixpar(rho = list(30, 30), eta = list(30, 30), xi = list(0.08,0.08),
##'              phi = list(c(K2 = 0.1, K3 = 0.2, K1 = 0.7), c(K2 = 0.1, K3 = 0.2, K1 = 0.7)))
##' rpmK <- rPeakHeight(mixK, pK, nsim = 1000, markers = "D2S1338")
##' oldpar <- par(mfrow = c(1,2))
##' boxplot(mixK, rpmK, markers = "D2S1338")
##' par(oldpar)
##' }
boxplot.DNAmixture <- function(x, sims, traces = 1:x$ntraces,
                               markers = x$markers, pw = 0.4,
                               ylim = NULL, border = "grey", ...){

  if (is.null(ylim)){
    ## Largest peak height
    ylim <- c(0,max(sapply(markers, function(m)max(x$data[[m]][,traces+1], sims[[m]][traces,,], na.rm = TRUE))))
  }

  bp <- function(m, r){
    d <- x$data[[m]]
    a <- d$allele

    heights <- sims[[m]][r,,]
    # colors <- apply(heights, 1, function(x)ifelse(all(x == 0), "blue", "grey"))
    pl <- plot(x, markers = m, traces = r, ylim = ylim, ask = FALSE, ...)
    allelepos <- pl[[m]][[r]]$mid

    usr <- par("usr")
    for (i in seq_along(allelepos)){
      clip(usr[1], usr[2], x$C[[r]], usr[4])
      boxplot(t(heights)[,i], at = allelepos[i], add = TRUE, border = border)
    }
    ## reset to plot region
    do.call("clip", as.list(usr))

    ## Re-draw peaks on top
    plot(x, markers = m, traces = r, add = TRUE, ask = FALSE,...)

    ## And points on top too
    points(allelepos, d[, r+1],
           col = "red", pch = ifelse(d[,r+1] > 0, 4, 1), cex = 1.3)
  }

  for (r in traces) {
    for (m in markers){
      if (r %in% x$observed[[m]]) bp(m, r)
    }
  }

  invisible(sims)
}
